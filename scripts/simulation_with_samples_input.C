#include <iostream>
#include <string>
#include <algorithm>
#include <tuple>
#include <vector>
#include <iterator>
#include <random>
#include <TRandom3.h>
#include <TTree.h>
#include <TFile.h>


/*****************************************************************************************************************
 * This script add to a truth level DAFNE sample the detector effect and background with the.
 * help of external samples (probably GEANT based samples). A truth level DAFNE generated samples
 * and two external samples are needed. The two external samples are:
 * 1. signal MC sample: this sample should include the signal process events. Each event should
 *                      contain both the truth level information and reconstruction level 
 *                      information. Then one can get the difference between the truth level 
 *                      and reconstructoin level for some variables (\delta v). With this \delta v,
 *                      one can produce a new reconstruction level variables with:
 *												             v_{rec} = v_{DAFNE} + \delta v
 *										  , where the v_{DAFNE} is the truth level variables from a events in the 
 *										  DAFNE sample. In this way, one simulated the resolution for variable v, nearly
 *										  as good as the GEANT simulation, without any modeling. For example, one don't
 *										  need to assume the resolution to be a gaussian or double gaussian distribution.
 *
 *											For the reconstruction efficiency, now we still need to model it into a weighting
 *											function. In the future development, one may feed a phase space GEANT based signal 
 *											MC sample into the DAFNE framework. Then the DAFNE framework may sample based on
 *											this reconstruction efficiency considered phase space. We can also combine the
 *											v_{DAFNE} and \delta v during this sampling procedure. In this way, we simulate both
 *											the efficiency and resolution without modeling.
 *
 * 2. generic MC sample: this sampel should include the background process events or signal and background
 *											 process events. Each event should contain the reconstruction level information
 *											 and the tag that indicating itself to be a signal or background event. Then one
 *											 can bootstrap the background events from this sample, without modeling the shape
 *											 of the background.
 *
 * The input samples could be large. This script can split them into small batches, while not losing randomness.
 *****************************************************************************************************************/



#define ISSIGNALTYPE bool 
// #define ISSIGNALTYPE double 

#define NLIMITVECTOR 5000000 // limit the size of the vector, to meet with the memory requirement of the queue
#define NROUND 10 

std::clock_t last_time_stamp = 0;

void TimeElapsedFor(std::string tag)
{
	if (tag == "") {
		std::cout << "start recording ... ..." << std::endl;
		last_time_stamp = std::clock(); 
		return;
	}
	
	std::clock_t current_time_stamp = std::clock();
	std::cout << "Time elapsed for " << tag << ": " << double(current_time_stamp - last_time_stamp) / CLOCKS_PER_SEC << " second" << std::endl;
	last_time_stamp = current_time_stamp;
	return;
}

struct Event {
	// explicit Event(double _m2p, double _m2m, double _t, double _sigmat, double _flavor, double _M, double _Q) : 
	//					m2p(_m2p), m2m(_m2m), t(_t), sigmat(_sigmat), flavor(_flavor), M(_M), Q(_Q) {};
	// Event() = default;

	// void print() const { }
	double m2p;
	double m2m;
	double t;
	double sigmat;
	double flavor;
	double M;
	double Q;
	double t_truth;
};

//Parameters:
//	dafneFilePath: the path of an ntuple file generated with dafne, with efficiency
//	genericMCFilePath: the path of an ntuple file contains candidates from genericMC
//	signalMCFilePath: the path of an ntuple file contains candidates from signalMC
void simulation_with_samples_input(const char * dafneFilePath, 
																const char * signalMCFilePath,
																const char * genericMCFilePath, 
																const char * outputFilePath,
																const double t_min, const double t_max, const double sigmat_min, const double sigmat_max, 
																int definedSignalYieldsInSignalRegion=-1, 
																int seed=10000, int signalOnly=0)
{


	std::string dafneTreeName = "ntp";
	std::string dafneM2pStr = "mSqP";
	std::string dafneM2mStr = "mSqM";
	std::string dafneM2zStr = "mSqZ";
	std::string dafneCosPiPiStr = "HelicityZ";
	std::string dafneDecayTimeStr = "t";
	std::string dafneSigmatStr = "sigmat";
	std::string dafneFlavorStr = "flavor";

	std::string genericMCTreeName = "Dst_KS_Nis";
	std::string genericMCMD0Str = "D0_M_nod0masscon";
	std::string genericMCQStr = "Q";
	std::string genericMCM2pStr = "m2p";
	std::string genericMCM2mStr = "m2m";
	std::string genericMCM2zStr = "m2z";
	std::string genericMCCosPiPiStr = "cospipi";
	std::string genericMCDecayTimeStr = "t";
	std::string signalMCTruthDecayTimeStr = "D0_mcFlightTime";
	std::string genericMCSigmatStr = "sigma_t";
	std::string genericMCFlavorStr = "spi_charge";
	std::string genericMCDst_isSignalStr = "Dst_isSignalAcceptMissingGammaAndDecayInFlight";
	std::string genericMCD0_isSignalStr = "D0_isSignalAcceptMissingGammaAndDecayInFlight";

	// **** input TTree preparation *****
	TFile * dafneFile = TFile::Open(dafneFilePath);
	TTree * dafneTree = (TTree*)dafneFile->Get(dafneTreeName.c_str()); 
	std::cout << "Finish loading dafne tree ......" << std::endl;
	TFile * genericMCFile = TFile::Open(genericMCFilePath);
	TTree * genericMCTree = (TTree*)genericMCFile->Get(genericMCTreeName.c_str()); 
	std::cout << "Finish loading generic MC tree ......" << std::endl;
	TFile * signalMCFile = TFile::Open(signalMCFilePath);
	TTree * signalMCTree = (TTree*)signalMCFile->Get(genericMCTreeName.c_str()); 
	std::cout << "Finish loading signal MC tree ......" << std::endl;

	// **** temp vectors *****
	std::vector<Event> dafneVector;
	std::vector<Event> signalMC_signalVector;
	std::vector<Event> signalMC_randomVector;
	std::vector<Event> signalVector;
	std::vector<Event> randomVector;
	std::vector<Event> otherVector;
		
	// **** output TTree preparation *****
	TFile * outputFile = TFile::Open(outputFilePath, "recreate");
	TTree * outputTree = new TTree(genericMCTreeName.c_str(), genericMCTreeName.c_str());

	const double M_low = 1.8;
	const double M_high = 2.0;
	const double Q_low = 0;
	const double Q_high = 20;

	const double MSR_low = 1.8648 - 0.015;
	const double MSR_high = 1.8648 + 0.015;
	const double QSR_low = 4.85;
	const double QSR_high = 6.85;

	
	double m2p_dafne;
	double m2m_dafne;
	double t_dafne;
	int flavor_dafne;

	double m2p_genericMC;
	double m2m_genericMC;
	double t_genericMC;
	double sigmat_genericMC;
	double flavor_genericMC;
	double M_genericMC;
	double Q_genericMC;
	ISSIGNALTYPE Dst_isSignal_genericMC;
	double D0_isSignal_genericMC;

	double t_signalMC;
	double t_truth_signalMC;
	double sigmat_signalMC;
	double flavor_signalMC;
	double M_signalMC;
	double Q_signalMC;
	ISSIGNALTYPE Dst_isSignal_signalMC;
	double D0_isSignal_signalMC;

	dafneTree->SetBranchAddress(dafneM2pStr.c_str(), &m2p_dafne);
	dafneTree->SetBranchAddress(dafneM2mStr.c_str(), &m2m_dafne);
	dafneTree->SetBranchAddress(dafneDecayTimeStr.c_str(), &t_dafne);
	dafneTree->SetBranchAddress(dafneFlavorStr.c_str(), &flavor_dafne);

	genericMCTree->SetBranchAddress(genericMCM2pStr.c_str(), &m2p_genericMC);
	genericMCTree->SetBranchAddress(genericMCM2mStr.c_str(), &m2m_genericMC);
	genericMCTree->SetBranchAddress(genericMCDecayTimeStr.c_str(), &t_genericMC);
	genericMCTree->SetBranchAddress(genericMCSigmatStr.c_str(), &sigmat_genericMC);
	genericMCTree->SetBranchAddress(genericMCFlavorStr.c_str(), &flavor_genericMC);
	genericMCTree->SetBranchAddress(genericMCMD0Str.c_str(), &M_genericMC);
	genericMCTree->SetBranchAddress(genericMCQStr.c_str(), &Q_genericMC);
	genericMCTree->SetBranchAddress(genericMCDst_isSignalStr.c_str(), &Dst_isSignal_genericMC);
	genericMCTree->SetBranchAddress(genericMCD0_isSignalStr.c_str(), &D0_isSignal_genericMC);

	signalMCTree->SetBranchAddress(genericMCDecayTimeStr.c_str(), &t_signalMC);
	signalMCTree->SetBranchAddress(signalMCTruthDecayTimeStr.c_str(), &t_truth_signalMC);
	signalMCTree->SetBranchAddress(genericMCSigmatStr.c_str(), &sigmat_signalMC);
	signalMCTree->SetBranchAddress(genericMCFlavorStr.c_str(), &flavor_signalMC);
	signalMCTree->SetBranchAddress(genericMCMD0Str.c_str(), &M_signalMC);
	signalMCTree->SetBranchAddress(genericMCQStr.c_str(), &Q_signalMC);
	signalMCTree->SetBranchAddress(genericMCDst_isSignalStr.c_str(), &Dst_isSignal_signalMC);
	signalMCTree->SetBranchAddress(genericMCD0_isSignalStr.c_str(), &D0_isSignal_signalMC);

	// read the data to temp vector, to speed up the bootstrapping
	TRandom3 r(seed);
	std::mt19937 rshuffle(seed);
	std::vector<Event> tempVector; // load the entries from TTree to tempVector sequentially. Sequential loading is faster then 
																 //	load with random index.
																 //	Then shuffle the tempVector and get the first NLIMITVECTOR entries, to save memory. Why we need 
																 //	shuffling? Because in some samples, the D0 and D0bar are not randomly distribute by indices. For example,
																 // the first half is all D0 and the other hald is all D0bar.
	double dafneTreeMaxTemp = dafneTree->GetEntries();
	// split to NROUND rounds, due to memory limitation
	double nmax_single_round = dafneTreeMaxTemp/NROUND;
	for (int i_round = 0; i_round < NROUND; ++i_round) {
		std::cout << "round " << i_round << std::endl;
		for (int i = i_round*nmax_single_round; i < (i_round+1)*nmax_single_round; ++i) {
			dafneTree->GetEntry(i);
			Event current_event;
			current_event.m2p = m2p_dafne;
			current_event.m2m = m2m_dafne;
			current_event.t = t_dafne*1000;
			current_event.flavor = flavor_dafne;
			tempVector.push_back(current_event);

			if (i%100000 == 0) std::cout << "reading dafneTree to temp vector" << i << std::endl;
		}

		std::cout << "Shuffling and loading into dafne vector ......" << std::endl;
		std::shuffle(tempVector.begin(), tempVector.end(), rshuffle);
		dafneVector.insert(dafneVector.end(), tempVector.begin(), tempVector.begin()+(NLIMITVECTOR/NROUND));
		tempVector.clear();
		tempVector.shrink_to_fit();
		std::cout << "Size of the dafne vector is " << dafneVector.size();

	}
	std::cout << "Finish loading dafne vector ......" << std::endl;


	for (int i = 0; i < genericMCTree->GetEntries(); ++i) {
		genericMCTree->GetEntry(i);
		Event current_event;
		current_event.m2p = m2p_genericMC;
		current_event.m2m = m2m_genericMC;
		current_event.t = t_genericMC;
		current_event.sigmat = sigmat_genericMC;
		current_event.flavor = flavor_genericMC;
		current_event.M = M_genericMC;
		current_event.Q = Q_genericMC;
		if (!(M_low<current_event.M&&current_event.M<M_high && Q_low<current_event.Q&&current_event.Q<Q_high)) continue;
		if (Dst_isSignal_genericMC == 1 && D0_isSignal_genericMC == 1)
			signalVector.push_back(current_event);	
		if (!(Dst_isSignal_genericMC == 1) && D0_isSignal_genericMC == 1) // !(isSignal==1) includes the nan case
			if (!signalOnly) randomVector.push_back(current_event);	
		if (!(Dst_isSignal_genericMC == 1) && !(D0_isSignal_genericMC == 1))
			if (!signalOnly) otherVector.push_back(current_event);	

		if (i%100000 == 0) std::cout << "reading genericMCTree " << i << std::endl;
	}

	std::cout << "Finish loading generic MC vector ......" << std::endl;

	std::vector<Event> temp_signalVector; 
	std::vector<Event> temp_randomVector;
	double signalMCTreeMaxTemp = signalMCTree->GetEntries();
	// split to NROUND rounds, due to memory limitation
	nmax_single_round = signalMCTreeMaxTemp/NROUND;
	for (int i_round = 0; i_round < NROUND; ++i_round) {
		std::cout << "round " << i_round << std::endl;
		for (int i = i_round*nmax_single_round; i < (i_round+1)*nmax_single_round; ++i) {
			signalMCTree->GetEntry(i);
			Event current_event;
			current_event.t = t_signalMC;
			current_event.t_truth = t_truth_signalMC*1000000;
			current_event.sigmat = sigmat_signalMC;
			current_event.flavor = flavor_signalMC;
			current_event.M = M_signalMC;
			current_event.Q = Q_signalMC;
			if (!(M_low<current_event.M&&current_event.M<M_high && Q_low<current_event.Q&&current_event.Q<Q_high)) continue;
			if (Dst_isSignal_signalMC == 1 && D0_isSignal_signalMC == 1)
				temp_signalVector.push_back(current_event);	
			if (!(Dst_isSignal_signalMC == 1) && D0_isSignal_signalMC == 1) // !(isSignal==1) includes the nan case
				temp_randomVector.push_back(current_event);	

			if (i%100000 == 0) std::cout << "reading signalMCTree to temp vectors" << i << std::endl;
		}

		std::cout << "Shuffling and loading into signal MC vectors ......" << std::endl;
		std::shuffle(temp_signalVector.begin(), temp_signalVector.end(), rshuffle);
		std::shuffle(temp_randomVector.begin(), temp_randomVector.end(), rshuffle);
		signalMC_signalVector.insert(signalMC_signalVector.end(), temp_signalVector.begin(), temp_signalVector.begin()+(NLIMITVECTOR/NROUND*(temp_signalVector.size())/(temp_signalVector.size()+temp_randomVector.size())));
		signalMC_randomVector.insert(signalMC_randomVector.end(), temp_randomVector.begin(), temp_randomVector.begin()+(NLIMITVECTOR/NROUND*(temp_randomVector.size())/(temp_signalVector.size()+temp_randomVector.size())));
		temp_signalVector.clear();
		temp_signalVector.shrink_to_fit();
		temp_randomVector.clear();
		temp_randomVector.shrink_to_fit();
		std::cout << "Size of the signalMC signal vector is " << signalMC_signalVector.size();
		std::cout << "Size of the signalMC random vector is " << signalMC_randomVector.size();
	}

	std::cout << "Finish loading signal MC vector ......" << std::endl;

	double n_signal = signalVector.size();
	double n_random = randomVector.size();
	double n_other = otherVector.size();
	double omega = 0.5; // mistag rate
	// double n_signal = 95000;
	// double n_random = 2000;
	// double n_other = 3000;
	double n_all = n_signal + n_random + n_other;

	double p_signal = n_signal/n_all;
	double p_random = (n_signal+n_random)/n_all;
	double p_other = (n_signal+n_random+n_other)/n_all;

	int n_gen = 0;
	if (definedSignalYieldsInSignalRegion != -1) {
		// assume the purity to be 95%
		double targetYieldsInSignalRegion = double(definedSignalYieldsInSignalRegion) / 0.95; 
		if (signalOnly) targetYieldsInSignalRegion = definedSignalYieldsInSignalRegion;
		int currentYieldsInSignalRegion = 0;
		for (int i_signal = 0; i_signal < signalVector.size(); ++i_signal) {
			Event & signal_event = signalVector[i_signal];
			if (!(MSR_low<signal_event.M&&signal_event.M<MSR_high)) continue;
			if (!(QSR_low<signal_event.Q&&signal_event.Q<QSR_high)) continue;
			currentYieldsInSignalRegion = currentYieldsInSignalRegion + 1;
		}
		std::cout << "ratio(signal yields in SR / all signal yields) = " << double(currentYieldsInSignalRegion)/double(n_signal) << " = " << int(currentYieldsInSignalRegion) << " / " << int(n_signal) << std::endl;
		double scaleFactor = double(targetYieldsInSignalRegion)/double(currentYieldsInSignalRegion);
		std::cout << "targetYieldsInSignalRegion = " << targetYieldsInSignalRegion << std::endl;
		std::cout << "scaleFactor = " << scaleFactor << std::endl;

		n_all = n_all*scaleFactor;
	}
	n_all = r.Poisson(n_all);

	printf("n_all = %d\n", int(n_all));
	printf("p_signal = %f\n", p_signal);
	printf("p_random = %f\n", p_random-p_signal);
	printf("p_other = %f\n", p_other-p_random);

	int nmax_dafne = dafneVector.size();
	int nmax_signalMC_signal = signalMC_signalVector.size();
	int nmax_signalMC_random = signalMC_randomVector.size();
	int nmax_signal = signalVector.size();
	int nmax_random = randomVector.size();
	int nmax_other = otherVector.size();

	double m2p_output;
	double m2m_output;
	double t_output;
	double sigmat_output;
	double flavor_output;
	double M_output;
	double Q_output;
	bool Dst_isSignal_output;
	bool D0_isSignal_output;
	bool issig_output;
	bool isrnd_output;
	bool isoth_output;

	outputTree->Branch(genericMCM2pStr.c_str(), &m2p_output, (genericMCM2pStr+"/D").c_str());
	outputTree->Branch(genericMCM2mStr.c_str(), &m2m_output, (genericMCM2mStr+"/D").c_str());
	outputTree->Branch(genericMCDecayTimeStr.c_str(), &t_output, (genericMCDecayTimeStr+"/D").c_str());
	outputTree->Branch(genericMCSigmatStr.c_str(), &sigmat_output, (genericMCSigmatStr+"/D").c_str());
	outputTree->Branch(genericMCFlavorStr.c_str(), &flavor_output, (genericMCFlavorStr+"/D").c_str());
	outputTree->Branch(genericMCMD0Str.c_str(), &M_output, (genericMCMD0Str+"/D").c_str());
	outputTree->Branch(genericMCQStr.c_str(), &Q_output, (genericMCQStr+"/D").c_str());
	outputTree->Branch(genericMCDst_isSignalStr.c_str(), &Dst_isSignal_output, (genericMCDst_isSignalStr+"/O").c_str());
	outputTree->Branch(genericMCD0_isSignalStr.c_str(), &D0_isSignal_output, (genericMCD0_isSignalStr+"/O").c_str());
	outputTree->Branch("issig", &issig_output, "issig/O");
	outputTree->Branch("isrnd", &isrnd_output, "isrnd/O");
	outputTree->Branch("isoth", &isoth_output, "isoth/O");


	while (n_gen < n_all) {
		// the current event is a signal candidate, random slow pion background candidate 
		// or other background candidate
		// TimeElapsedFor("");
		double p_category = r.Uniform(1);
		int category = 0; // 1 for signal, 2 for random slow background, 3 for other background
		if (p_category < p_signal) category = 1;
		else if (p_category < p_random) category = 2;
		else category = 3;
		// TimeElapsedFor("category decision");

		if (category == 1) {
			int i_signal = r.Integer(nmax_signal);
			Event & signal_event = signalVector[i_signal];

			Event & dafne_event = dafneVector[r.Integer(nmax_dafne)];	
			flavor_output = dafne_event.flavor;

			Event & signalMC_signal_event = signalMC_signalVector[r.Integer(nmax_signalMC_signal)];
			while (signalMC_signal_event.flavor != flavor_output)	signalMC_signal_event = signalMC_signalVector[r.Integer(nmax_signalMC_signal)];

			// TimeElapsedFor("category 1: reading entries from trees");
			m2p_output = dafne_event.m2p;
			m2m_output = dafne_event.m2m;

			// add some weighting code for the reconstruction efficiency here
			// ................weighting for the reconstruction efficiency...............
			// Or, when generating the truth level samples with the DAFNE, you can add also the efficiency simulation,
			// then you don't need to weight here.


			t_output = dafne_event.t + (signalMC_signal_event.t - signalMC_signal_event.t_truth);
			sigmat_output = signalMC_signal_event.sigmat;
			M_output = signalMC_signal_event.M;
			Q_output = signalMC_signal_event.Q;

			Dst_isSignal_output = 1;
			D0_isSignal_output = 1;
			
			issig_output = 1;
			isrnd_output = 0;
			isoth_output = 0;

			// TimeElapsedFor("category 1: filling tree");
		} else if (category == 2) {
			// **** random slow pion background ******
			// ** bootstrapping Q distribution
			int i_random = r.Integer(nmax_random);
			Event & random_event = randomVector[i_random];
			Event & signalMC_random_event = signalMC_randomVector[r.Integer(nmax_signalMC_random)];
			Event & dafne_event = dafneVector[r.Integer(nmax_dafne)];	
			// TimeElapsedFor("category 2: reading entries from trees");
			m2p_output = dafne_event.m2p;
			m2m_output = dafne_event.m2m;
			t_output = dafne_event.t + (signalMC_random_event.t - signalMC_random_event.t_truth);
			sigmat_output = signalMC_random_event.sigmat;
			flavor_output = dafne_event.flavor;

			// One may change the mistag rate omega here
			// ........... updating the omega ................

			if (r.Uniform(1) < omega) {
				flavor_output = -dafne_event.flavor;
				m2p_output = dafne_event.m2m;
				m2m_output = dafne_event.m2p;
			}
			M_output = random_event.M; // currently, for the soft pion, we directly use the M and Q from the generic MC,
			Q_output = random_event.Q; // ignoring the correlation between MQ and decay-time variables.
			
			Dst_isSignal_output = 0;
			D0_isSignal_output = 1;

			issig_output = 0;
			isrnd_output = 1;
			isoth_output = 0;
			
			// TimeElapsedFor("category 2: filling tree");
		}	else if (category == 3) {
			// **** other background ******
			// ** bootstrapping
			int i_other = r.Integer(nmax_other);
			// TimeElapsedFor("category 3: reading entries from trees");
			Event & other_event = otherVector[i_other];
			m2p_output = other_event.m2p;
			m2m_output = other_event.m2m;
			t_output = other_event.t;
			sigmat_output = other_event.sigmat;
			flavor_output = other_event.flavor;
			M_output = other_event.M;
			Q_output = other_event.Q;

			Dst_isSignal_output = 0;
			D0_isSignal_output = 0;

			issig_output = 0;
			isrnd_output = 0;
			isoth_output = 1;

			// TimeElapsedFor("category 3: filling tree");
		} else {
			std::cout << "wrong category : " << category << std::endl;
			exit(0);
		}

		if (!(t_min<t_output && t_output<t_max && sigmat_min<sigmat_output && sigmat_output<sigmat_max)) continue;
		outputTree->Fill();

		
		n_gen++;
		// if (n_gen%100 == 0) std::cout << n_gen << std::endl;
		
	} 


	outputFile->cd();
	outputTree->Write();
	outputFile->Close();
	
}

