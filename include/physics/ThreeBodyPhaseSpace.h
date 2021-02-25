#pragma once

#include <random>

#include <hydra/Lambda.h>
#include <hydra/PhaseSpace.h>
#include <hydra/Decays.h>
#include <hydra/SparseHistogram.h>
#include <hydra/SeedRNG.h>
#include <hydra/Vector4R.h>
#include <hydra/functions/Exponential.h>
#include <hydra/functions/UniformShape.h>
#include <hydra/functions/Gaussian.h>
#include <hydra/Plain.h>
#include <hydra/GenzMalikQuadrature.h>


#include <TH3D.h>
#include <THnSparse.h>

#include <physics/Utils.h>
#include <physics/Rate.h>

namespace dafne {

class ThreeBodyPhaseSpace
{
private:
	double _mMother;
	std::array<double,3> _mDaughters;
	
public:
	ThreeBodyPhaseSpace() = delete;

	__hydra_dual__
	ThreeBodyPhaseSpace(double const &mass_mother, const std::array<double,3> &mass_daugthers) : _mMother(mass_mother), _mDaughters(mass_daugthers) {}

	__hydra_dual__
	ThreeBodyPhaseSpace(double const &mass_mother, double const &mass_daugther1, double const &mass_daugther2, double const &mass_daugther3) : ThreeBodyPhaseSpace(mass_mother,{mass_daugther1,mass_daugther2,mass_daugther3}) {}
   
	__hydra_dual__ inline
	ThreeBodyPhaseSpace(ThreeBodyPhaseSpace const& other) : _mMother(other.MMother()), _mDaughters(other.MDaughters()) {}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpace& operator=(ThreeBodyPhaseSpace const& other)
	{
		if(this == &other) return *this;
		_mMother = other.MMother();
		_mDaughters = other.MDaughters();
		return *this;
	}
	
	__hydra_dual__ inline
	double MMother() const { return _mMother; }

	__hydra_dual__ inline
   double MSqMother() const { return _mMother*_mMother; }

	__hydra_dual__ inline
	std::array<double,3> MDaughters() const { return _mDaughters; }
	 
	template<unsigned Index>
	__hydra_dual__ inline
	double M() const
	{
		HYDRA_STATIC_ASSERT(Index<=3 && Index>=1,"Index must be in [1,3].");
		return _mDaughters[Index-1];
	}
	
	template<unsigned Index>
	__hydra_dual__ inline
	double MSq() const { return M<Index>()*M<Index>(); }
	
	__hydra_dual__ inline
	double MSqSum() const { return MSqMother()+MSq<1>()+MSq<2>()+MSq<3>(); }

	__hydra_dual__ inline
	double MSqJK(double const &mSqIJ, double const &mSqIK) const { return MSqSum() - mSqIJ - mSqIK; };
	

	// cos(theta(Helicity)) of the JK system
	template<unsigned J, unsigned K>
	__hydra_dual__ inline
	double HelicityJK(double const &mSqIJ, double const &mSqIK) const 
	{
		double mSqJK = MSqJK(mSqIJ, mSqIK);
		double mJK = sqrt(mSqJK);
		double erK = (mSqJK - std::pow(M<J>(),2) + std::pow(M<K>(),2)) / (2*mJK); // energy in the rest frame
		double erI = (std::pow(_mMother,2) - mSqJK - std::pow(M<6-J-K>(),2)) / (2*mJK);
		double prK = sqrt(std::pow(erK,2) - std::pow(M<K>(),2)); // momentum in the rest frame
		double prI = sqrt(std::pow(erI,2) - std::pow(M<6-J-K>(),2));
		double smaxIK = std::pow((erK + erI),2) - std::pow((prK - prI),2);
		double sminIK = std::pow((erK + erI),2) - std::pow((prK + prI),2);
		double helicityJK = (smaxIK + sminIK - 2*mSqIK) / (smaxIK - sminIK);
		return helicityJK;
	}

	// convert HelicityJK back to mSqIJ, mSqIK
	template<unsigned J, unsigned K>
	__hydra_dual__ inline
	auto invHelicityJK(double const &helicityJK, double const &mSqJK) const 
	{
		double edrJK = (std::pow(_mMother,2) - std::pow(M<6-J-K>(),2) + mSqJK)/(2*_mMother); // energy in D0 rest frame
		double mJK = sqrt(mSqJK);
		double pdrJK = sqrt(std::pow(edrJK,2) - std::pow(mJK,2)); // momentum in D0 rest frame
		double edrI = _mMother - edrJK;
		double pdrI = sqrt(std::pow(edrI,2) - std::pow(M<6-J-K>(),2));
		double beta = pdrJK / edrJK;
		double gamma = 1 / sqrt(1 - beta*beta);
		double erI = gamma*(edrI + beta*pdrI);
		double prI = sqrt(std::pow(erI,2) - std::pow(M<6-J-K>(),2)); // momentum in JK rest frame
		double erK = mJK / 2;
		double prK = sqrt(std::pow(erK,2) - std::pow(M<K>(),2));
		double erJ = erK;
		double prJ = prK;
		double mSqIK = std::pow(M<6-J-K>(),2) + std::pow(M<K>(),2) + 2*erI*erK - 2*prI*prK*helicityJK;
		double mSqIJ = std::pow(_mMother,2) + std::pow(M<6-J-K>(),2) + std::pow(M<J>(),2) + std::pow(M<K>(),2) - mSqIK - mSqJK;

		return std::make_pair(mSqIJ, mSqIK);
	}

	template<unsigned I, unsigned J>
	__hydra_dual__ inline
	double MSqMin() const { return std::pow( M<I>() + M<J>(), 2 ); };
	
	template<unsigned I, unsigned J>
	__hydra_dual__ inline
	double MSqMax() const { return std::pow( _mMother - M<6-I-J>(), 2 ); };
	
	template<unsigned I, unsigned J>
	__hydra_dual__ inline
	double MSqMin( const double& mSqIK ) const
	{
		double first  = MSqMother() + MSq<I>() - MSq<J>() - MSq<6-I-J>();
		double second = std::sqrt( kallen( mSqIK, MSq<I>()   , MSq<6-I-J>() ) );
		second       += std::sqrt( kallen( mSqIK, MSqMother(), MSq<J>()     ) );

		return ( first*first - second*second ) / ( 4. * mSqIK );
	}

	template<unsigned I, unsigned J>
	__hydra_dual__ inline
	double MSqMax( const double& mSqIK ) const
	{
		double first  = MSqMother() + MSq<I>() - MSq<J>() - MSq<6-I-J>();
		double second = std::sqrt( kallen( mSqIK, MSq<I>()   , MSq<6-I-J>() ) );
		second       -= std::sqrt( kallen( mSqIK, MSqMother(), MSq<J>()     ) );

		return ( first*first - second*second ) / ( 4. * mSqIK );
	}
	
	template<unsigned J, unsigned K>
	__hydra_dual__ inline
	double MSqMin( const double& mSqIJ, const double& mSqIK ) const
	{
		double firstIJ  = MSqMother() - MSq<6-J-K>() + MSq<J>() - MSq<K>();
		double secondIJ = std::sqrt( kallen( mSqIJ, MSq<6-J-K>(), MSq<J>() ) );
		secondIJ       += std::sqrt( kallen( mSqIJ, MSqMother() , MSq<K>() ) );

		double firstIK  = MSqMother() - MSq<6-J-K>() - MSq<J>() + MSq<K>();
		double secondIK = std::sqrt( kallen( mSqIK, MSq<6-J-K>(), MSq<K>() ) );
		secondIK       += std::sqrt( kallen( mSqIK, MSqMother() , MSq<J>() ) );

		return std::max( ( firstIJ*firstIJ - secondIJ*secondIJ ) / ( 4. * mSqIJ ),
						     ( firstIK*firstIK - secondIK*secondIK ) / ( 4. * mSqIK ) );
	}
	
	template<unsigned J, unsigned K>
	__hydra_dual__ inline
	double MSqMax( const double& mSqIJ, const double& mSqIK ) const
	{
		double firstIJ  = MSqMother() - MSq<6-J-K>() + MSq<J>() - MSq<K>();
		double secondIJ = std::sqrt( kallen( mSqIJ, MSq<6-J-K>(), MSq<J>() ) );
		secondIJ       -= std::sqrt( kallen( mSqIJ, MSqMother() , MSq<K>() ) );

		double firstIK  = MSqMother() - MSq<6-J-K>() - MSq<J>() + MSq<K>();
		double secondIK = std::sqrt( kallen( mSqIK, MSq<6-J-K>(), MSq<K>() ) );
		secondIK       -= std::sqrt( kallen( mSqIK, MSqMother() , MSq<J>() ) );

		return std::max( ( firstIJ*firstIJ - secondIJ*secondIJ ) / ( 4. * mSqIJ ),
					   	  ( firstIK*firstIK - secondIK*secondIK ) / ( 4. * mSqIK ) );
	}
	
	template<unsigned J, unsigned K>
	__hydra_dual__ inline
	bool Contains( const double& mSqIJ, const double& mSqIK ) const
	{
		const double& mSqJK = MSqSum() - mSqIJ - mSqIK;
		return ( mSqIK > MSqMin<6-J-K,K>( mSqIJ ) ) && ( mSqIK < MSqMax<6-J-K,K>( mSqIJ ) ) &&
				 ( mSqJK > MSqMin<J,K>    ( mSqIJ ) ) && ( mSqJK < MSqMax<J,K>    ( mSqIJ ) );
	}

	__hydra_dual__ inline
	bool Contains( const double& mSq12, const double& mSq13, const double& mSq23 ) const
	{
		return ( mSq12 > MSqMin<1,2>( mSq13        ) ) && ( mSq12 < MSqMax<1,2>( mSq13        ) ) &&
				 ( mSq13 > MSqMin<1,3>( mSq12        ) ) && ( mSq13 < MSqMax<1,3>( mSq12        ) ) &&
				 ( mSq23 > MSqMin<2,3>( mSq12, mSq13 ) ) && ( mSq23 < MSqMax<2,3>( mSq12, mSq13 ) );
	}
	
	template<typename RND=hydra::default_random_engine>
	__hydra_dual__ inline
	auto Generator() { return hydra::PhaseSpace<3,RND>{_mMother,_mDaughters}; }
	
	template<typename A=hydra::Vector4R, typename B=hydra::Vector4R, typename C=hydra::Vector4R, typename Backend=hydra::device::sys_t>
	__hydra_dual__ inline
	auto Decays(size_t n) { return hydra::Decays<hydra::tuple<A,B,C>, Backend>(_mMother,_mDaughters,n); }

	template<typename RND=hydra::default_random_engine>
	__hydra_dual__ inline
	auto Integrator(size_t nevents=1000000)
	{
		nevents = (nevents>=1000000) ? nevents : 1000000; // force nevents = 1000000 if nevents is too small

		std::array<double,2> min{MSqMin<1,2>(), MSqMin<1,3>()};
		std::array<double,2> max{MSqMax<1,2>(), MSqMax<1,3>()};
		return hydra::Plain<2, hydra::device::sys_t, RND>(min,max,nevents);
	}
	
	__hydra_dual__ inline
	THnSparseD *RootHistogram(const char *name1="1", const char *name2="2", const char *name3="3", int nbins=300)
	{
		const int bins[] = {nbins, nbins, nbins};
		const double min[] = {MSqMin<1,2>(),MSqMin<1,3>(),MSqMin<2,3>()};
		const double max[] = {MSqMax<1,2>(),MSqMax<1,3>(),MSqMax<2,3>()};
		auto histo = new THnSparseD("histo","", 3, bins, min, max);
		histo->GetAxis(0)->SetTitle(Form("#it{m}^{2}(%s%s) [GeV^{2}/#it{c}^{4}]",name1,name2));
		histo->GetAxis(1)->SetTitle(Form("#it{m}^{2}(%s%s) [GeV^{2}/#it{c}^{4}]",name1,name3));
		histo->GetAxis(2)->SetTitle(Form("#it{m}^{2}(%s%s) [GeV^{2}/#it{c}^{4}]",name2,name3));
		return histo;
	}
	
	template<typename Backend=hydra::device::sys_t>
	__hydra_dual__ inline
	auto SparseHistogram(size_t nbins=300)
	{
		return hydra::SparseHistogram<double, 3, Backend>{
			{nbins, nbins, nbins},{MSqMin<1,2>(),MSqMin<1,3>(),MSqMin<2,3>()},{MSqMax<1,2>(),MSqMax<1,3>(),MSqMax<2,3>()}
		};
	}
	
	template<typename MSq12, typename MSq13, typename MSq23, typename Model>
	__hydra_dual__ inline
	auto GenerateData(Model &model, size_t nevents, size_t rndseed=0)
	{
		// Output container
		hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23>, hydra::device::sys_t> data;
		
		// Functor to compute Dalitz variables from 4-momenta
		auto dalitz_calculator = hydra::wrap_lambda( [] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			MSq23 m2 = (b+c).mass2();
			return hydra::make_tuple(m0,m1,m2);
		});
		
		// Model weighting functor
		auto dalitz_model = hydra::wrap_lambda( [model] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			return model(hydra::tie(m0,m1));
		});
		
		// Create PhaseSpace object
		auto phsp_generator = Generator();
		
		// Sample total number of events to be generated
		hydra::SeedRNG seed{rndseed};
		
		static std::mt19937_64 random_mt(seed());
		std::poisson_distribution<size_t> poisson(nevents);
		size_t n = poisson(random_mt);
		
		// Generate in bunches
		auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(10*n);
		hydra::Vector4R parent(MMother(), 0., 0., 0.);
		do {
			phsp_generator.SetSeed(seed());
			phsp_generator.Generate(parent, phsp_events);
				
			// Unweight
			auto dalitz_variables = phsp_events.Unweight(dalitz_model, -1, seed()) | dalitz_calculator;
				
			// First copy unweighted events into a container
			hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23>, hydra::device::sys_t> bunch( dalitz_variables.size() );
			hydra::copy(dalitz_variables, bunch);
				
			// Then add to output container
			data.insert( data.end(), bunch.begin(), bunch.end() );
		} while( data.size() < n );
		
		// Erase excess of events
		data.erase(data.begin()+n, data.end());
		return data;
	}

	template<typename MSq12, typename MSq13, typename MSq23, typename Model>
	__hydra_dual__ inline
	auto GenerateSparseHistogram(Model &model, size_t nevents, size_t nbins=300, size_t rndseed=0)
	{
		// Functor to compute Dalitz variables from 4-momenta
		auto dalitz_calculator = hydra::wrap_lambda( [] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			MSq23 m2 = (b+c).mass2();
			return hydra::make_tuple(m0,m1,m2);
		});
		
		// Model weighting functor
		auto dalitz_model = hydra::wrap_lambda( [model] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			return model(hydra::tie(m0,m1));
		});
		
		// Generate flat phase-space
		hydra::SeedRNG seed{rndseed};
		hydra::Vector4R parent(MMother(), 0., 0., 0.);
		auto phsp_generator = Generator();
		phsp_generator.SetSeed(seed());
		auto events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(nevents);
		phsp_generator.Generate(parent, events);
		
		// Compute Dalitz variables and weights
		auto variables = events | dalitz_calculator;
		auto weights   = events | events.GetEventWeightFunctor(dalitz_model);
		
		// Fill histogram
		auto histo = SparseHistogram(nbins);
		histo.Fill(variables, weights);
		
		return histo;
	}
	

	// Compute fit fractions
	template<typename MSq12, typename MSq13, typename Model, typename Amplitude>
	__hydra_dual__ inline
	auto FitFraction(Model const& model, Amplitude const& amp, size_t nentries=1000000)
	{
		auto integrator = Integrator(nentries);
		
		auto ampsq = hydra::wrap_lambda( [amp] __hydra_dual__ (MSq12 a, MSq13 b) {
			return hydra::norm( amp(hydra::tie(a,b)) );
		});

		auto modelsq = hydra::wrap_lambda( [model] __hydra_dual__ (MSq12 a, MSq13 b) {
			return hydra::norm( model(hydra::tie(a,b)) );
		});

		hydra::SeedRNG seed{};
		integrator.SetSeed(seed());
		auto int_acomp = integrator.Integrate(ampsq);
		integrator.SetSeed(seed()); 
		auto int_model = integrator.Integrate(modelsq);
	
		auto f = int_acomp.first/int_model.first;
		auto e = f * sqrt( pow(int_acomp.second/int_acomp.first,2) + pow(int_model.second/int_model.first,2) );
		return std::make_pair(f,e);
	}

};

	
class ThreeBodyPhaseSpaceWithTime : public ThreeBodyPhaseSpace
{
private:
	std::array<double,2> _timeRange;
	
public:
	ThreeBodyPhaseSpaceWithTime() = delete;
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTime(ThreeBodyPhaseSpace const &phsp, std::array<double,2> const &time_range) : ThreeBodyPhaseSpace(phsp), _timeRange(time_range) {}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTime(double const &mass_mother, const std::array<double,3> &mass_daugthers, std::array<double,2> const &time_range) : ThreeBodyPhaseSpace(mass_mother,mass_daugthers), _timeRange(time_range) {}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTime(double const &mass_mother, double const &mass_daugther1, double const &mass_daugther2, double const &mass_daugther3, double const &t_min, double const &t_max) : ThreeBodyPhaseSpace(mass_mother,{mass_daugther1,mass_daugther2,mass_daugther3}), _timeRange({t_min,t_max}) {}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTime(ThreeBodyPhaseSpaceWithTime const& other) : ThreeBodyPhaseSpaceWithTime(other.MMother(),other.MDaughters(),other.TimeRange()) {}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTime& operator=(ThreeBodyPhaseSpaceWithTime const& other)
	{
		if(this == &other) return *this;
		ThreeBodyPhaseSpace::operator=(other);
		_timeRange = other.TimeRange();
		return *this;
	}
	
	__hydra_dual__ inline
	std::array<double,2> TimeRange() const { return _timeRange; }

	__hydra_dual__ inline
	double TimeMin() const { return _timeRange[0]; }
	
	__hydra_dual__ inline
	double TimeMax() const { return _timeRange[1]; }
	
	template<typename RND=hydra::default_random_engine>
	__hydra_dual__ inline
	auto IntegratorWithTime(size_t nevents=1000000)
	{
		nevents = (nevents>=1000000) ? nevents : 1000000; // force nevents = 1000000 if nevents is too small
		
		std::array<double,3> min{TimeMin(), MSqMin<1,2>(), MSqMin<1,3>()};
		std::array<double,3> max{TimeMax(), MSqMax<1,2>(), MSqMax<1,3>()};
		return hydra::Plain<3, hydra::device::sys_t, RND>(min,max,nevents);
	}
	
	__hydra_dual__ inline
	THnSparseD *RootHistogramWithTime(const char *name1="1", const char *name2="2", const char *name3="3", int nbins=300)
	{
		const int bins[] = {nbins, nbins, nbins, nbins};
		const double min[] = {MSqMin<1,2>(),MSqMin<1,3>(),MSqMin<2,3>(),TimeMin()};
		const double max[] = {MSqMax<1,2>(),MSqMax<1,3>(),MSqMax<2,3>(),TimeMax()};
		auto histo = new THnSparseD("histo","", 4, bins, min, max);
		histo->GetAxis(0)->SetTitle(Form("#it{m}^{2}(%s%s) [GeV^{2}/#it{c}^{4}]",name1,name2));
		histo->GetAxis(1)->SetTitle(Form("#it{m}^{2}(%s%s) [GeV^{2}/#it{c}^{4}]",name1,name3));
		histo->GetAxis(2)->SetTitle(Form("#it{m}^{2}(%s%s) [GeV^{2}/#it{c}^{4}]",name2,name3));
		histo->GetAxis(3)->SetTitle("#it{t} [ps]");
		return histo;
	}

	template<typename Backend=hydra::device::sys_t>
	__hydra_dual__ inline
	auto SparseHistogramWithTime(size_t nbins=300)
	{
		return hydra::SparseHistogram<double, 4, Backend>{
			{nbins, nbins, nbins, nbins},{MSqMin<1,2>(),MSqMin<1,3>(),MSqMin<2,3>(),TimeMin()},{MSqMax<1,2>(),MSqMax<1,3>(),MSqMax<2,3>(),TimeMax()}
		};
	}

	template<typename MSq12, typename MSq13, typename MSq23, typename Time, typename Model>
	__hydra_dual__ inline
	auto GenerateDataWithTime(Model &model, size_t nevents, size_t rndseed=0, double y=-999, double Gamma=-999)
	{
		// Output container
		hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23,Time>, hydra::device::sys_t> data;

		// Functor to compute Dalitz variables from 4-momenta
		auto dalitz_calculator = hydra::wrap_lambda( [] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			MSq23 m2 = (b+c).mass2();
			return hydra::make_tuple(m0,m1,m2,t);
		});

		// Create PhaseSpace object
		auto phsp_generator = Generator();
		hydra::Vector4R parent(MMother(), 0., 0., 0.);
		
		// Find maximum
		double max_model(-1.);
		{
			auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(10*nevents);
			phsp_generator.Generate(parent, phsp_events);
			auto phsp_weight = phsp_events.GetEventWeightFunctor();
		
			auto dalitz_model = hydra::wrap_lambda( [&phsp_weight, &model] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c)
			{
				MSq12 m0 = (a+b).mass2();
				MSq13 m1 = (a+c).mass2();
				Time t(0.);
				return phsp_weight(a,b,c) * model(hydra::tie(t,m0,m1));
			});
		
			auto variables = phsp_events | dalitz_model;
			max_model = *( hydra_thrust::max_element(hydra::device::sys, variables.begin(), variables.end() ) );
		}

		// Sample total number of events to be generated
		hydra::SeedRNG seed{rndseed};

		static std::mt19937_64 random_mt(seed());
		std::poisson_distribution<size_t> poisson(nevents);
		size_t n = poisson(random_mt);
				
		// Generated in bunches
		auto uniform = hydra::UniformShape<Time>(TimeMin(),TimeMax());
		auto uniform_for_exp_outline = hydra::UniformShape<double>(0,1);
		// the inverse function of cdf of pdf(t) = (1-|y|) Gamma * e^[ -(1-|y|) Gamma t ], 
		// used only when !(y!=-999 && Gamma!=-999), to sample time according to Jordi's 
		// method
		auto exp_outline_invcdf = hydra::wrap_lambda( [y, Gamma] __hydra_dual__ (double q)
		{
			Time t = std::log(1-q) / ( -(1.0-abs(y))*Gamma )  ;
			return t;
		});

		hydra::device::vector<Time> time_data(10*n);
		hydra::device::vector<Time> uniform_data(time_data.size());
		auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(time_data.size());
		do {
			phsp_generator.SetSeed(seed());
			phsp_generator.Generate(parent, phsp_events);

			// Phase-space weighting functor
			auto phsp_weight = phsp_events.GetEventWeightFunctor();

			// Model weighting functor
			auto dalitz_time_model = hydra::wrap_lambda( [&phsp_weight, &model, y, Gamma] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t)
			{
				MSq12 m0 = (a+b).mass2();
				MSq13 m1 = (a+c).mass2();
				double weight = phsp_weight(a,b,c) * model(hydra::tie(m0,m1,t));
				if (y!=-999 && Gamma!=-999) return weight * 1. / ((1.-abs(y))*Gamma * exp(-(1-abs(y))*Gamma*t)); // weight under exponential outline
				return weight;
			});

			// generated decay time uniformuly when Gamma and y is not set, otherwise generate with the exp(-(1-abs(y))*Gamma*t)
			// just like Jordi's code
			if (!(y!=-999 && Gamma!=-999)) hydra::fill_random( time_data, uniform, seed()); 
			else {
				hydra::fill_random( uniform_data, uniform_for_exp_outline, seed()); // the RngFormula<Exponential> looks confusing, 
				hydra::copy(uniform_data | exp_outline_invcdf, time_data);	    	// I would like to use the explicit invcdf lambda above
			}
			

			auto events_with_time = phsp_events.Meld( time_data );

			// Unweight
			auto dalitz_variables = hydra::unweight(hydra::device::sys, events_with_time, dalitz_time_model, 1.1*max_model, seed()) | dalitz_calculator;

			// First copy unweighted events into a container
			hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23,Time>, hydra::device::sys_t> bunch( dalitz_variables.size() );
			hydra::copy(dalitz_variables, bunch);

			// Then add to output container
			data.insert( data.end(), bunch.begin(), bunch.end() );
		} while( data.size() < n );

		// Erase excess of events
		data.erase(data.begin()+n, data.end());
		return data;
	}

	template<typename MSq12, typename MSq13, typename MSq23, typename Time, typename Model>
	__hydra_dual__ inline
	auto GenerateSparseHistogramWithTime(Model &model, size_t nevents, size_t nbins=300, size_t rndseed=0)
	{
		// Functor to compute Dalitz variables from 4-momenta
		auto dalitz_calculator = hydra::wrap_lambda( [] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			MSq23 m2 = (b+c).mass2();
			return hydra::make_tuple(m0,m1,m2,t);
		});

		// Generate flat phase-space
		hydra::SeedRNG seed{rndseed};
		hydra::Vector4R parent(MMother(), 0., 0., 0.);
		auto phsp_generator = Generator();
		phsp_generator.SetSeed(seed());
		auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(nevents);
		phsp_generator.Generate(parent, phsp_events);

		// Add uniformuly generated decay time to phase-space data
		hydra::device::vector<Time> time_data(nevents);
		hydra::fill_random( time_data, hydra::UniformShape<Time>(TimeMin(),TimeMax()), seed());

		auto events = phsp_events.Meld( time_data );

		// Phase-space weighting functor
		auto phsp_weight = phsp_events.GetEventWeightFunctor();

		// Model weighting functor
		auto dalitz_time_model = hydra::wrap_lambda( [phsp_weight, model] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			return phsp_weight(a,b,c) * model(hydra::tie(t,m0,m1));
		});

		// Compute Dalitz variables and weights
		auto variables = events | dalitz_calculator;
		auto weights   = events | dalitz_time_model;

		// Fill histogram
		auto histo = SparseHistogramWithTime(nbins);
		histo.Fill(variables, weights);

		return histo;
	}
};


class ThreeBodyPhaseSpaceWithTimeAndTimeError : public ThreeBodyPhaseSpace
{
private:
	std::array<double,2> _timeRange;
	std::array<double,2> _timeErrorRange;
	
public:
	ThreeBodyPhaseSpaceWithTimeAndTimeError() = delete;
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTimeAndTimeError(ThreeBodyPhaseSpace const &phsp, std::array<double,2> const &time_range, std::array<double,2> const &time_error_range) : ThreeBodyPhaseSpace(phsp), _timeRange(time_range), _timeErrorRange(time_error_range) {}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTimeAndTimeError(double const &mass_mother, const std::array<double,3> &mass_daugthers, std::array<double,2> const &time_range, std::array<double,2> const &time_error_range) : ThreeBodyPhaseSpace(mass_mother,mass_daugthers), _timeRange(time_range), _timeErrorRange(time_error_range) {}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTimeAndTimeError(double const &mass_mother, double const &mass_daugther1, double const &mass_daugther2, double const &mass_daugther3, double const &t_min, double const &t_max, double const &sigmat_min, double const &sigmat_max) : ThreeBodyPhaseSpace(mass_mother,{mass_daugther1,mass_daugther2,mass_daugther3}), _timeRange({t_min,t_max}), _timeErrorRange({sigmat_min, sigmat_max}) {}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTimeAndTimeError(ThreeBodyPhaseSpaceWithTimeAndTimeError const& other) : ThreeBodyPhaseSpaceWithTimeAndTimeError(other.MMother(),other.MDaughters(),other.TimeRange(),other.TimeErrorRange()) {}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpaceWithTimeAndTimeError& operator=(ThreeBodyPhaseSpaceWithTimeAndTimeError const& other)
	{
		if(this == &other) return *this;
		ThreeBodyPhaseSpace::operator=(other);
		_timeRange = other.TimeRange();
		_timeErrorRange = other.TimeErrorRange();
		return *this;
	}
	
	__hydra_dual__ inline
	std::array<double,2> TimeRange() const { return _timeRange; }

	__hydra_dual__ inline
	std::array<double,2> TimeErrorRange() const { return _timeErrorRange; }

	__hydra_dual__ inline
	double TimeMin() const { return _timeRange[0]; }
	
	__hydra_dual__ inline
	double TimeMax() const { return _timeRange[1]; }

	__hydra_dual__ inline
	double TimeErrorMin() const { return _timeErrorRange[0]; }
	
	__hydra_dual__ inline
	double TimeErrorMax() const { return _timeErrorRange[1]; }

	__hydra_dual__ inline
	auto GenzMalikIntegratorWithTimeAndTimeError(size_t grid=30)
	// still has segment fault problem, and needs further check
	{
		std::array<double,4> min{MSqMin<1,2>(), MSqMin<1,3>(), TimeMin(), TimeErrorMin()};
		std::array<double,4> max{MSqMax<1,2>(), MSqMax<1,3>(), TimeMax(), TimeErrorMax()};
		return hydra::GenzMalikQuadrature<4, hydra::device::sys_t>(min,max,grid);
	}

	template<typename RND=hydra::default_random_engine>
	__hydra_dual__ inline
	auto PlainIntegratorWithTimeAndTimeError(size_t nevents=1000000)
	{
		nevents = (nevents>=1000000) ? nevents : 1000000; // force nevents = 1000000 if nevents is too small
		
		std::array<double,4> min{MSqMin<1,2>(), MSqMin<1,3>(), TimeMin(), TimeErrorMin()};
		std::array<double,4> max{MSqMax<1,2>(), MSqMax<1,3>(), TimeMax(), TimeErrorMax()};
		return hydra::Plain<4, hydra::device::sys_t, RND>(min,max,nevents);
	}
	
	__hydra_dual__ inline
	THnSparseD *RootHistogramWithTimeAndTimeError(const char *name1="1", const char *name2="2", const char *name3="3", int nbins=300)
	{
		const int bins[] = {nbins, nbins, nbins, nbins, nbins};
		const double min[] = {MSqMin<1,2>(),MSqMin<1,3>(),MSqMin<2,3>(),TimeMin(),TimeErrorMin()};
		const double max[] = {MSqMax<1,2>(),MSqMax<1,3>(),MSqMax<2,3>(),TimeMax(),TimeErrorMax()};
		auto histo = new THnSparseD("histo","", 5, bins, min, max);
		histo->GetAxis(0)->SetTitle(Form("#it{m}^{2}(%s%s) [GeV^{2}/#it{c}^{4}]",name1,name2));
		histo->GetAxis(1)->SetTitle(Form("#it{m}^{2}(%s%s) [GeV^{2}/#it{c}^{4}]",name1,name3));
		histo->GetAxis(2)->SetTitle(Form("#it{m}^{2}(%s%s) [GeV^{2}/#it{c}^{4}]",name2,name3));
		histo->GetAxis(3)->SetTitle("#it{t} [ps]");
		histo->GetAxis(4)->SetTitle("#it{#sigma_{t}} [ps]");
		return histo;
	}

	template<typename Backend=hydra::device::sys_t>
	__hydra_dual__ inline
	auto SparseHistogramWithTimeAndTimeError(size_t nbins=300)
	{
		return hydra::SparseHistogram<double, 5, Backend>{
			{nbins, nbins, nbins, nbins, nbins},{MSqMin<1,2>(),MSqMin<1,3>(),MSqMin<2,3>(),TimeMin(),TimeErrorMin()},{MSqMax<1,2>(),MSqMax<1,3>(),MSqMax<2,3>(),TimeMax(),TimeErrorMax()}
		};
	}

	template<typename MSq12, typename MSq13, typename MSq23, typename Time, typename TimeError, typename Model>
	__hydra_dual__ inline
	auto GenerateDataWithTimeAndTimeError(Model &model, size_t nevents, size_t rndseed=0, double TimeMaxSearchLowBound=-999., double TimeMaxSearchHighBound=-999., double TimeErrorMaxSearchLowBound=-999., double TimeErrorMaxSearchHighBound=-999., bool debug=false)
	{
		// Output container
		hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23,Time,TimeError>, hydra::device::sys_t> data;

		// Functor to compute Dalitz variables from 4-momenta
		auto dalitz_calculator = hydra::wrap_lambda( [] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t, TimeError sigmat)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			MSq23 m2 = (b+c).mass2();
			return hydra::make_tuple(m0,m1,m2,t,sigmat);
		});

		// Create PhaseSpace object
		auto phsp_generator = Generator();
		hydra::Vector4R parent(MMother(), 0., 0., 0.);
		
		// initialize the random seed
		hydra::SeedRNG seed{rndseed};


		// Find maximum

		// For the moment, the maximum value of sigmat distribution (which is JohnsonSU function in many cases) 
		// is currently searched by numerical method, and the range of the maximum value search could be set by input 
		// arguments.
		// The sampling time performance could be improved, if analytical method is used.

		// After considering the resolution, model(t = 0) is not longer max(model(t)), therefore, the maximum value 
		// of decay-time distribution should also be searched by numerical method, and the range of the search could
		// also be set by input arguments.
		double max_model(-1.);
		{
			auto max_search_seed = seed();
			int max_search_nevents = (nevents*5 > 10000000) ? nevents*5 : 10000000;
			auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(max_search_nevents);
			hydra::multivector< hydra::tuple<Time, TimeError>, hydra::device::sys_t > time_data(max_search_nevents);

			phsp_generator.Generate(parent, phsp_events);
			auto phsp_weight = phsp_events.GetEventWeightFunctor();

			if (TimeMaxSearchLowBound != -999. && TimeMaxSearchHighBound != -999.) 
				hydra::fill_random(time_data.begin(hydra::placeholders::_0), time_data.end(hydra::placeholders::_0), hydra::UniformShape<Time>(TimeMaxSearchLowBound, TimeMaxSearchHighBound), seed());	
			else
				hydra::fill_random(time_data.begin(hydra::placeholders::_0), time_data.end(hydra::placeholders::_0), hydra::UniformShape<Time>(TimeMin(), TimeMax()), seed());

			if (TimeErrorMaxSearchLowBound != -999. && TimeErrorMaxSearchHighBound != -999.) 
				hydra::fill_random(time_data.begin(hydra::placeholders::_1), time_data.end(hydra::placeholders::_1), hydra::UniformShape<TimeError>(TimeErrorMaxSearchLowBound, TimeErrorMaxSearchHighBound), max_search_seed);
			else 
				hydra::fill_random(time_data.begin(hydra::placeholders::_1), time_data.end(hydra::placeholders::_1), hydra::UniformShape<TimeError>(TimeErrorMin(), TimeErrorMax()), max_search_seed);

			

			auto events = phsp_events.Meld( time_data );

			auto dalitz_model = hydra::wrap_lambda( [&phsp_weight, &model] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t, TimeError sigmat)
			{
				MSq12 m0 = (a+b).mass2();
				MSq13 m1 = (a+c).mass2();
				return phsp_weight(a,b,c) * model(hydra::tie(m0,m1,t,sigmat));
			});
		
			auto variables = events | dalitz_model;
			max_model = *( hydra_thrust::max_element(hydra::device::sys, variables.begin(), variables.end() ) );

			if (debug) {
				std::cout << "Finish the search for max value of model(mSqP, mSqM, t, sigma_t). " << std::endl; 
				std::cout << "max(model(mSqP, mSqM, t, sigma_t)) = " << max_model << std::endl;
			}

		}

		// Sample total number of events to be generated
		static std::mt19937_64 random_mt(seed());
		std::poisson_distribution<size_t> poisson(nevents);
		size_t n = poisson(random_mt);

		// Generated in bunches
		long bunch_size = 10*n;
		long bunch_size_limit = 50000000;
		if (300*n < bunch_size_limit) bunch_size = 300*n;
		else if (100*n < bunch_size_limit) bunch_size = 100*n;
		else if (50*n < bunch_size_limit) bunch_size = 50*n;
		else if (10*n < bunch_size_limit) bunch_size = 10*n;
		else bunch_size = 3*n;

		hydra::multivector< hydra::tuple<Time, TimeError>, hydra::device::sys_t > time_data(bunch_size);

		auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(time_data.size());
		do {
			phsp_generator.SetSeed(seed());
			phsp_generator.Generate(parent, phsp_events);

			// Phase-space weighting functor
			auto phsp_weight = phsp_events.GetEventWeightFunctor();

			// Model weighting functor
			auto dalitz_time_model = hydra::wrap_lambda( [&phsp_weight, &model] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t, TimeError sigmat)
			{
				MSq12 m0 = (a+b).mass2();
				MSq13 m1 = (a+c).mass2();
				return phsp_weight(a,b,c) * model(hydra::tie(m0,m1,t,sigmat));
			});

			// Add uniformuly generated decay time and decay time error to phase-space data
			hydra::fill_random(time_data.begin(hydra::placeholders::_0), time_data.end(hydra::placeholders::_0), hydra::UniformShape<Time>(TimeMin(), TimeMax()), seed());
			hydra::fill_random(time_data.begin(hydra::placeholders::_1), time_data.end(hydra::placeholders::_1), hydra::UniformShape<TimeError>(TimeErrorMin(), TimeErrorMax()), seed());

			auto events = phsp_events.Meld( time_data );

			// confirm there is no pdfmax value issue
			if (debug) {
				double bunch_max_model(-1.);
				auto pdfmax_check_variables = events | dalitz_time_model;
				bunch_max_model = *( hydra_thrust::max_element(hydra::device::sys, pdfmax_check_variables.begin(), pdfmax_check_variables.end() ) );
				std::cout << "Events number for the current bunch = " << pdfmax_check_variables.size() << std::endl;
				std::cout << "For the current bunch, max(model(mSqP, mSqM, t, sigma_t)) = " << bunch_max_model << std::endl; 
			}

			// Unweight
			auto dalitz_variables = hydra::unweight(hydra::device::sys, events, dalitz_time_model, 1.3*max_model, seed()) | dalitz_calculator;

			// First copy unweighted events into a container
			hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23,Time,TimeError>, hydra::device::sys_t> bunch( dalitz_variables.size() );
			hydra::copy(dalitz_variables, bunch);

			if (debug) std::cout << "Generated " << bunch.size() << " events" << std::endl; 

			// Then add to output container
			data.insert( data.end(), bunch.begin(), bunch.end() );
		} while( data.size() < n );

		// Erase excess of events
		data.erase(data.begin()+n, data.end());
		return data;
	}

	// This function first generate (mSqP, mSqM, t) in the same way as GenerateDataWithTime(), then it sample a sigmat from pdf(sigmat) and smear
	// t with Gauss(t, b, s*sigmat).
	// Therefore it is a 3D + 1D sampling, which would be much faster than the 4D sampling in GenerateDataWithTimeAndTimeError and easier
	// to find a proper pdfmax.
	template<typename MSq12, typename MSq13, typename MSq23, typename Time, typename TimeError, typename ModelTruth, typename PDFSIGMAT>
	__hydra_dual__ inline
	auto GenerateDataWithTimeAndTimeErrorFast1(ModelTruth &modelTruth, double b, double s, PDFSIGMAT const& pdf_sigma_t,  size_t nevents, size_t rndseed=0,  double y=-999, double Gamma=-999, bool debug=false)
	{
			// Output container
		hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23,Time,TimeError>, hydra::device::sys_t> data;


		// Create PhaseSpace object
		auto phsp_generator = Generator();
		hydra::Vector4R parent(MMother(), 0., 0., 0.);
		
		// Find maximum
		double max_model(-1.);
		{
			auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(10*nevents);
			phsp_generator.Generate(parent, phsp_events);
			auto phsp_weight = phsp_events.GetEventWeightFunctor();
		
			auto dalitz_model = hydra::wrap_lambda( [&phsp_weight, &modelTruth] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c)
			{
				MSq12 m0 = (a+b).mass2();
				MSq13 m1 = (a+c).mass2();
				Time t(0.);
				return phsp_weight(a,b,c) * modelTruth(hydra::tie(t,m0,m1));
			});
		
			auto variables = phsp_events | dalitz_model;
			max_model = *( hydra_thrust::max_element(hydra::device::sys, variables.begin(), variables.end() ) );


			if (debug) {
				std::cout << "Finish the search for max value of model(mSqP, mSqM, t=0). " << std::endl; 
				std::cout << "max(model(mSqP, mSqM, t=0)) = " << max_model << std::endl;
			}
		}

		// Sample total number of events to be generated
		hydra::SeedRNG seed{rndseed};

		static std::mt19937_64 random_mt(seed());
		std::poisson_distribution<size_t> poisson(nevents);
		size_t n = poisson(random_mt);
				
		// Generated in bunches
		auto uniform = hydra::UniformShape<Time>(0,TimeMax()*10); // truth level decay-time axis should start from 0
		auto uniform_for_exp_outline = hydra::UniformShape<double>(0,1); // information of this sampling method could be found in GenerateDataWithTime()
		auto exp_outline_invcdf = hydra::wrap_lambda( [y, Gamma] __hydra_dual__ (double q)
		{
			Time t = std::log(1-q) / ( -(1.0-abs(y))*Gamma )  ;
			return t;
		});



		auto gaus = hydra::Gaussian<double>(b,s);

		hydra::multivector< hydra::tuple<Time, TimeError, double> , hydra::device::sys_t> time_data(10*n);
		hydra::device::vector<Time> uniform_data(time_data.size());
		auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(time_data.size());


		// Functor to compute Dalitz variables from 4-momenta
		auto dalitz_calculator = hydra::wrap_lambda( [] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t, TimeError sigma_t, double dtDsigmat)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			MSq23 m2 = (b+c).mass2();
			TimeError tsmear = dtDsigmat * sigma_t;
			Time trec = t + tsmear;
			return hydra::make_tuple(m0,m1,m2,trec,sigma_t);
		});

		do {
			phsp_generator.SetSeed(seed());
			phsp_generator.Generate(parent, phsp_events);

			// Phase-space weighting functor
			auto phsp_weight = phsp_events.GetEventWeightFunctor();

			// Model weighting functor
			auto dalitz_time_model = hydra::wrap_lambda( [&phsp_weight, &modelTruth, y, Gamma, this] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t, TimeError sigma_t, double dtDsigmat)
			{
				// judge in range
				TimeError tsmear = dtDsigmat * sigma_t;
				Time trec = t + tsmear;
				if (!(TimeMin()<trec && trec<TimeMax())) return 0.0;
				if (!(TimeErrorMin()<sigma_t && sigma_t<TimeErrorMax())) return 0.0;

				MSq12 m0 = (a+b).mass2();
				MSq13 m1 = (a+c).mass2();
				double weight = phsp_weight(a,b,c) * modelTruth(hydra::tie(m0,m1,t));

				if (y!=-999 && Gamma!=-999) return weight * 1. / ((1.-abs(y))*Gamma * exp(-(1-abs(y))*Gamma*t)); // weight under exponential outline
				return weight;
			});

			// Add uniformuly generated decay time to phase-space data

			if (!(y!=-999 && Gamma!=-999)) hydra::fill_random(time_data.begin(hydra::placeholders::_0), time_data.end(hydra::placeholders::_0), uniform, seed());
			else {
				hydra::fill_random( uniform_data, uniform_for_exp_outline, seed()); 
				auto time_data_temp = uniform_data | exp_outline_invcdf;
				hydra::copy(time_data_temp.begin(), time_data_temp.end(), time_data.begin(hydra::placeholders::_0));	    	
			}
			

			// sample sigma_t 
			hydra::fill_random(time_data.begin(hydra::placeholders::_1), time_data.end(hydra::placeholders::_1), pdf_sigma_t, seed());

			// sample deltat/sigma_t
			hydra::fill_random(time_data.begin(hydra::placeholders::_2), time_data.end(hydra::placeholders::_2), gaus, seed());

			auto events_with_time = phsp_events.Meld( time_data );

			// confirm there is no pdfmax value issue
			if (debug) {
				double bunch_max_model(-1.);
				auto pdfmax_check_variables = events_with_time | dalitz_time_model;
				bunch_max_model = *( hydra_thrust::max_element(hydra::device::sys, pdfmax_check_variables.begin(), pdfmax_check_variables.end() ) );
				std::cout << "Events number for the current bunch = " << pdfmax_check_variables.size() << std::endl;
				std::cout << "For the current bunch, max(model(mSqP, mSqM, t)) = " << bunch_max_model << std::endl; 
			}

			// Unweight
			auto dalitz_variables = hydra::unweight(hydra::device::sys, events_with_time, dalitz_time_model, 1.1*max_model, seed()) | dalitz_calculator;

			// First copy unweighted events into a container
			hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23,Time,TimeError>, hydra::device::sys_t> bunch( dalitz_variables.size() );
			hydra::copy(dalitz_variables, bunch);

			// Then add to output container
			data.insert( data.end(), bunch.begin(), bunch.end() );
		} while( data.size() < n );

		// Erase excess of events
		data.erase(data.begin()+n, data.end());
		return data;
	}


	// This function first generate (mSqP, mSqM) in the same way as GenerateData(), then it sample t from 
	// a sigmat from pdf(sigmat) and smear
	// t with Gauss(t, b, s*sigmat).
	// Therefore it is a 3D + 1D sampling, which would be much faster than the 4D sampling in GenerateDataWithTimeAndTimeError and easier
	// to find a proper pdfmax.
	template<typename MSq12, typename MSq13, typename MSq23, typename Time, typename TimeError, typename ModelTruth, typename PDFSIGMAT>
	__hydra_dual__ inline
	auto GenerateDataWithTimeAndTimeErrorFast2(ModelTruth &modelTruth, double b, double s, PDFSIGMAT const& pdf_sigma_t,  size_t nevents, size_t rndseed=0, bool debug=false)
	{
			// Output container
		hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23,Time,TimeError>, hydra::device::sys_t> data;


		// Create PhaseSpace object
		auto phsp_generator = Generator();
		hydra::Vector4R parent(MMother(), 0., 0., 0.);
		
		// Find maximum
		double max_model(-1.);
		{
			auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(10*nevents);
			phsp_generator.Generate(parent, phsp_events);
			auto phsp_weight = phsp_events.GetEventWeightFunctor();
		
			auto dalitz_model = hydra::wrap_lambda( [&phsp_weight, &modelTruth] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c)
			{
				MSq12 m0 = (a+b).mass2();
				MSq13 m1 = (a+c).mass2();
				Time t(0.);
				return phsp_weight(a,b,c) * modelTruth(hydra::tie(t,m0,m1));
			});
		
			auto variables = phsp_events | dalitz_model;
			max_model = *( hydra_thrust::max_element(hydra::device::sys, variables.begin(), variables.end() ) );
		}

		// Sample total number of events to be generated
		hydra::SeedRNG seed{rndseed};

		static std::mt19937_64 random_mt(seed());
		std::poisson_distribution<size_t> poisson(nevents);
		size_t n = poisson(random_mt);
				
		// Generated in bunches
		auto uniform = hydra::UniformShape<Time>(0,TimeMax()*1.3); // truth level decay-time axis should start from 0
		auto gaus = hydra::Gaussian<double>(b,s);

		hydra::multivector< hydra::tuple<Time, TimeError, double> , hydra::device::sys_t> time_data(10*n);
		auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(time_data.size());


		// Functor to compute Dalitz variables from 4-momenta
		auto dalitz_calculator = hydra::wrap_lambda( [] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t, TimeError sigma_t, double dtDsigmat)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			MSq23 m2 = (b+c).mass2();
			TimeError tsmear = dtDsigmat * sigma_t;
			Time trec = t + tsmear;
			return hydra::make_tuple(m0,m1,m2,trec,sigma_t);
		});

		do {
			phsp_generator.SetSeed(seed());
			phsp_generator.Generate(parent, phsp_events);

			// Phase-space weighting functor
			auto phsp_weight = phsp_events.GetEventWeightFunctor();

			// Model weighting functor
			auto dalitz_time_model = hydra::wrap_lambda( [&phsp_weight, &modelTruth, this] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t, TimeError sigma_t, double dtDsigmat)
			{
				// judge in range
				TimeError tsmear = dtDsigmat * sigma_t;
				Time trec = t + tsmear;
				if (!(TimeMin()<t && t<TimeMax())) return 0.0;
				if (!(TimeErrorMin()<sigma_t && sigma_t<TimeErrorMin())) return 0.0;

				MSq12 m0 = (a+b).mass2();
				MSq13 m1 = (a+c).mass2();
				return phsp_weight(a,b,c) * modelTruth(hydra::tie(m0,m1,t));
			});

			// Add uniformuly generated decay time to phase-space data
			hydra::fill_random(time_data.begin(hydra::placeholders::_0), time_data.end(hydra::placeholders::_0), uniform, seed());

			// sample sigma_t 
			hydra::fill_random(time_data.begin(hydra::placeholders::_1), time_data.end(hydra::placeholders::_1), pdf_sigma_t, seed());

			// sample deltat/sigma_t
			hydra::fill_random(time_data.begin(hydra::placeholders::_2), time_data.end(hydra::placeholders::_2), gaus, seed());

			auto events_with_time = phsp_events.Meld( time_data );

			// Unweight
			auto dalitz_variables = hydra::unweight(hydra::device::sys, events_with_time, dalitz_time_model, 1.1*max_model, seed()) | dalitz_calculator;

			// First copy unweighted events into a container
			hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23,Time,TimeError>, hydra::device::sys_t> bunch( dalitz_variables.size() );
			hydra::copy(dalitz_variables, bunch);

			// Then add to output container
			data.insert( data.end(), bunch.begin(), bunch.end() );
		} while( data.size() < n );

		// Erase excess of events
		data.erase(data.begin()+n, data.end());
		return data;
	}


	template<typename MSq12, typename MSq13, typename MSq23, typename Time, typename TimeError, typename Model>
	__hydra_dual__ inline
	auto GenerateSparseHistogramWithTimeAndTimeError(Model &model, size_t nevents, size_t nbins=300, size_t rndseed=0)
	{
		// Functor to compute Dalitz variables from 4-momenta
		auto dalitz_calculator = hydra::wrap_lambda( [] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t, TimeError sigmat)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			MSq23 m2 = (b+c).mass2();
			return hydra::make_tuple(m0,m1,m2,t,sigmat);
		});

		// Generate flat phase-space
		hydra::SeedRNG seed{rndseed};
		hydra::Vector4R parent(MMother(), 0., 0., 0.);
		auto phsp_generator = Generator();
		phsp_generator.SetSeed(seed());
		auto phsp_events = Decays<hydra::Vector4R,hydra::Vector4R,hydra::Vector4R>(nevents);
		phsp_generator.Generate(parent, phsp_events);

		// Add uniformuly generated decay time to phase-space data
		hydra::multivector< hydra::tuple<Time, TimeError> , hydra::device::sys_t> time_data(nevents);
		hydra::fill_random(time_data.begin(hydra::placeholders::_0), time_data.end(hydra::placeholders::_0), hydra::UniformShape<Time>(TimeMin(),TimeMax()), seed());
		hydra::fill_random(time_data.begin(hydra::placeholders::_1), time_data.end(hydra::placeholders::_1), hydra::UniformShape<TimeError>(TimeErrorMin(),TimeErrorMax()), seed());

		auto events = phsp_events.Meld( time_data );

		// Phase-space weighting functor
		auto phsp_weight = phsp_events.GetEventWeightFunctor();

		// Model weighting functor
		auto dalitz_time_model = hydra::wrap_lambda( [phsp_weight, model] __hydra_dual__ (hydra::Vector4R a, hydra::Vector4R b, hydra::Vector4R c, Time t, TimeError sigmat)
		{
			MSq12 m0 = (a+b).mass2();
			MSq13 m1 = (a+c).mass2();
			return phsp_weight(a,b,c) * model(hydra::tie(m0,m1,t,sigmat));
		});

		// Compute Dalitz variables and weights
		auto variables = events | dalitz_calculator;
		auto weights   = events | dalitz_time_model;

		// Fill histogram
		auto histo = SparseHistogramWithTimeAndTimeError(nbins);
		histo.Fill(variables, weights);

		return histo;
	}
};




}//dafne namespace
