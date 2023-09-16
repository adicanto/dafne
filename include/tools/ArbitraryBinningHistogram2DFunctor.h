#ifndef ARBITRARYBINNINGHISTOGRAM2DFUNCTOR_H_
#define ARBITRARYBINNINGHISTOGRAM2DFUNCTOR_H_

#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Spiline.h>
#include <hydra/Placeholders.h>
#include <hydra/detail/utility/CheckValue.h>
#include <hydra/detail/external/hydra_thrust/copy.h>
#include <hydra/detail/external/hydra_thrust/iterator/zip_iterator.h>
#include <hydra/detail/external/hydra_thrust/execution_policy.h>
#include <hydra/detail/external/hydra_thrust/binary_search.h>
#include <hydra/detail/external/hydra_thrust/extrema.h>
#include <hydra/detail/external/hydra_thrust/iterator/iterator_traits.h>
#include <hydra/detail/external/hydra_thrust/functional.h>
#include <hydra/detail/external/hydra_thrust/advance.h>
#include <math.h>
#include <algorithm>
#include <memory>

#include <tools/ArbitraryBinningHistogram2D.h>

namespace hydra {


template<typename TypeZ, typename Iterator1, typename Iterator2, typename Iterator3, typename TypeX, typename TypeY>
__hydra_host__ __hydra_device__
inline typename std::enable_if< std::is_floating_point<typename hydra::thrust::iterator_traits<Iterator1>::value_type >::value &&
                       std::is_floating_point<typename hydra::thrust::iterator_traits<Iterator2>::value_type >::value && std::is_floating_point<typename hydra::thrust::iterator_traits<Iterator3>::value_type >::value, TypeZ>::type
histogram2d(Iterator1 firstXTick, Iterator1 lastXTick,  
		Iterator2 firstYTick, Iterator2 lastYTick,  
	    Iterator3 Z, TypeX valueX, TypeY valueY) {

		// std::cout << "search start   " << valueX <<","<< valueY << std::endl;
		auto iterX = hydra::detail::spiline::lower_bound(firstXTick, lastXTick, valueX);
		auto iterY = hydra::detail::spiline::lower_bound(firstYTick, lastYTick, valueY);
		size_t dist_i_x = hydra::thrust::distance(firstXTick, iterX);
		size_t dist_i_y = hydra::thrust::distance(firstYTick, iterY);
		size_t i_x = dist_i_x > 0 ? dist_i_x - 1: 0;
		size_t i_y = dist_i_y > 0 ? dist_i_y - 1: 0;

		// std::cout << "dist_i_x = " << dist_i_x << std::endl;
		// std::cout << "dist_i_y = " << dist_i_y << std::endl;
		// std::cout << "i_x = " << i_x << std::endl;
		// std::cout << "i_y = " << i_y << std::endl;

		size_t N_x = hydra::thrust::distance(firstXTick, lastXTick);
		size_t N_y = hydra::thrust::distance(firstYTick, lastYTick);

		// notice! the bin number is tick number - 1
		N_x = N_x-1;
		N_y = N_y-1;

		// std::cout << "N_x = " << N_x << std::endl;
		// std::cout << "N_y = " << N_y << std::endl;
		// std::cout << "i_y*N_x + i_x = " << i_y*N_x + i_x << std::endl;

		//--------------------

		const double z_i = Z[i_y*N_x + i_x];

		return z_i;
	}


template<typename TypeZ, typename Iterable1, typename Iterable2, typename Iterable3, typename TypeX, typename TypeY>
__hydra_host__ __hydra_device__
inline typename std::enable_if< hydra::detail::is_iterable<Iterable1>::value &&
                       hydra::detail::is_iterable<Iterable2>::value &&
                        hydra::detail::is_iterable<Iterable3>::value &&
                       std::is_floating_point<typename Iterable1::value_type >::value &&
                       std::is_floating_point<typename Iterable2::value_type >::value &&
                       std::is_floating_point<typename Iterable3::value_type >::value ,
                       TypeZ >::type
histogram2d(Iterable1&& XTicks, Iterable2&& YTicks, Iterable3&& Zs, TypeX valueX, TypeY valueY){

	return histogram2d<TypeZ>( std::forward<Iterable1>(XTicks).begin(),
			            std::forward<Iterable1>(XTicks).end(),
			            std::forward<Iterable2>(YTicks).begin(),
			            std::forward<Iterable2>(YTicks).end(),
						std::forward<Iterable3>(Zs).begin() , valueX, valueY);

}



}




namespace hydra {

template<typename Iterator1, typename Iterator2, typename Iterator3, typename Iterator4, typename ArgTypeX, typename ArgTypeY, typename ArgTypeZ, typename Signature=ArgTypeZ(ArgTypeX,ArgTypeY)>
class ArbitraryBinningHistogram2DFunctor: public BaseFunctor<ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, Iterator3, Iterator4, ArgTypeX, ArgTypeY, ArgTypeZ>, Signature, 0>
{


// A functor for arbitrary binning 2d histogram.
// fXTick: iterator for the edges of x bins
// fYTick: iterator for the edges of y bins
// fZ: iterator for the values in bins, the length of the bin value iterable should be
//     (fSizeXTick-1)*(fSizeYTick-1)
// fZError: iterator for the value errors in bins, the length of the bin value error iterable should be
//     (fSizeXTick-1)*(fSizeYTick-1)


public:

	ArbitraryBinningHistogram2DFunctor() = delete;

	ArbitraryBinningHistogram2DFunctor( Iterator1 firstXTick, Iterator1 lastXTick, 
										Iterator2 firstYTick, Iterator2 lastYTick, 
										Iterator3 Z, Iterator4 ZError 
										):
		BaseFunctor<ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, Iterator3, Iterator4, ArgTypeX, ArgTypeY, ArgTypeZ>, Signature, 0>(),
		fSizeXTick(hydra::thrust::distance(firstXTick, lastXTick)),
		fSizeYTick(hydra::thrust::distance(firstYTick, lastYTick)),
		fXTick(firstXTick),
		fYTick(firstYTick),
		fZ(Z),
		fZError(ZError)
		{}

	__hydra_host__ __hydra_device__
	ArbitraryBinningHistogram2DFunctor(ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, Iterator3, Iterator4, ArgTypeX, ArgTypeY, ArgTypeZ> const& other ):
	BaseFunctor<ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, Iterator3, Iterator4, ArgTypeX, ArgTypeY, ArgTypeZ>, Signature, 0>(other),
	fSizeXTick(other.GetSizeXTick()),
	fSizeYTick(other.GetSizeYTick()),
	fXTick(other.GetXTick()),
	fYTick(other.GetYTick()),
	fZ(other.GetZ()),
	fZError(other.GetZError())
	{ }

	__hydra_host__ __hydra_device__ inline
	ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, Iterator3, Iterator4, ArgTypeX, ArgTypeY, ArgTypeZ>&
	operator=(ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, Iterator3, Iterator4, ArgTypeX, ArgTypeY, ArgTypeZ> const& other )
	{
		if(this == &other) return *this;

		BaseFunctor<ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, Iterator3, Iterator4, ArgTypeX, ArgTypeY, ArgTypeZ>, Signature, 0>::operator=(other);

		fSizeXTick=other.GetSizeXTick();
		fSizeYTick=other.GetSizeYTick();
		fXTick=other.GetXTick();
		fYTick=other.GetYTick();
		fZ=other.GetZ();
		fZError=other.GetZError();

		return *this;
	}

	__hydra_host__ __hydra_device__
	size_t GetSizeXTick() const
	{
		return fSizeXTick;
	}

	__hydra_host__ __hydra_device__
	size_t GetSizeYTick() const
	{
		return fSizeYTick;
	}


	__hydra_host__ __hydra_device__
	Iterator1 GetXTick() const
	{
		return fXTick;
	}

	__hydra_host__ __hydra_device__
	Iterator2 GetYTick() const
	{
		return fYTick;
	}

	__hydra_host__ __hydra_device__
	Iterator3 GetZ() const
	{
		return fZ;
	}

	__hydra_host__ __hydra_device__
	Iterator4 GetZError() const
	{
		return fZError;
	}


	__hydra_host__ __hydra_device__
	inline double Evaluate(ArgTypeX x, ArgTypeY y)  const {


		Iterator1 fXTickN = fXTick + fSizeXTick;
		Iterator1 fYTickN = fYTick + fSizeYTick;

		double _x = x;
		double _y = y;

		// add range judegement here
		double z = histogram2d<double>(fXTick, fXTickN, fYTick, fYTickN, fZ, _x, _y);

		// currently only return the Z value
		// the returning of Z error is to be decided
		return  CHECK_VALUE( z, "z=%f", z) ;
	}



private:

	size_t fSizeXTick; 
	size_t fSizeYTick;

	Iterator1 fXTick;
	Iterator2 fYTick;

	Iterator3 fZ;
	Iterator4 fZError;

};


}  // namespace hydra



namespace hydra {

// template<typename ArgType, typename Iterator1, typename Iterator2, typename Iterator3, typename Iterator4>
// inline ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, ArgType>
// make_arbitrary_binning_histogram_2d_functor(Iterator1 firstXTick, Iterator1 lastXTick, 
// 											Iterator2 firstYTick, Iterator3 firstZ, 
// 											Iterator4 firstZError )
// {
// }

template<typename ArgTypeX, typename ArgTypeY, typename ArgTypeZ, typename Iterable1, typename Iterable2, typename Iterable3, typename Iterable4>
inline typename std::enable_if<
          hydra::detail::is_iterable<Iterable1>::value && hydra::detail::is_iterable<Iterable2>::value &&
          hydra::detail::is_iterable<Iterable3>::value && hydra::detail::is_iterable<Iterable4>::value ,
          ArbitraryBinningHistogram2DFunctor< decltype(std::declval<Iterable1>().begin()),
                          					  decltype(std::declval<Iterable2>().begin()), 
                          					  decltype(std::declval<Iterable3>().begin()), 
                          					  decltype(std::declval<Iterable4>().begin()), 
                          					  ArgTypeX, ArgTypeY, ArgTypeZ> >::type
make_arbitrary_binning_histogram_2d_functor(Iterable1&& XTicks, Iterable2&& YTicks, Iterable3&& Zs, Iterable4&& ZErrors)
{

typedef  decltype(std::declval<Iterable1>().begin()) Iterator1;
typedef  decltype(std::declval<Iterable2>().begin()) Iterator2;
typedef  decltype(std::declval<Iterable3>().begin()) Iterator3;
typedef  decltype(std::declval<Iterable4>().begin()) Iterator4;

	return ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, Iterator3, Iterator4, ArgTypeX, ArgTypeY, ArgTypeZ>(
			std::forward<Iterable1>(XTicks).begin(),
			std::forward<Iterable1>(XTicks).end(),
			std::forward<Iterable1>(YTicks).begin(),
			std::forward<Iterable1>(YTicks).end(),
			std::forward<Iterable2>(Zs).begin(),
			std::forward<Iterable2>(ZErrors).begin());
}

template<typename ArgTypeX, typename ArgTypeY, typename ArgTypeZ, typename Iterable1, typename Iterable2, typename Iterable3, typename Iterable4>
inline   typename std::enable_if< 
          hydra::detail::is_iterable<Iterable1>::value && hydra::detail::is_iterable<Iterable2>::value &&
          hydra::detail::is_iterable<Iterable3>::value && hydra::detail::is_iterable<Iterable4>::value ,
          ArbitraryBinningHistogram2DFunctor< decltype(std::declval<Iterable1>().begin()),
                          					  decltype(std::declval<Iterable2>().begin()), 
                          					  decltype(std::declval<Iterable3>().begin()), 
                          					  decltype(std::declval<Iterable4>().begin()), 
                          					 ArgTypeX, ArgTypeY, ArgTypeZ> >::type
make_arbitrary_binning_histogram_2d_functor(Iterable1 & XTicks, Iterable2 & YTicks, Iterable3 & Zs, Iterable4 & ZErrors, ArbitraryBinningHistogram2D & input_hist)
{

typedef  decltype(std::declval<Iterable1>().begin()) Iterator1;
typedef  decltype(std::declval<Iterable2>().begin()) Iterator2;
typedef  decltype(std::declval<Iterable3>().begin()) Iterator3;
typedef  decltype(std::declval<Iterable4>().begin()) Iterator4;

	XTicks.clear();
	YTicks.clear();
	Zs.clear();
	ZErrors.clear();

	std::cout << std::endl;
	for (int i_x = 0; i_x < input_hist.GetNXTicks(); ++i_x) {
		double v = input_hist.GetXTick(i_x);
		std::cout << v << ", "; 
		XTicks.push_back(v);
	}
	std::cout << std::endl;
	for (int i_x = 0; i_x < input_hist.GetNXTicks(); ++i_x) std::cout << XTicks[i_x] << ", "; 
	std::cout << std::endl;

	std::cout << std::endl;
	for (int i_y = 0; i_y < input_hist.GetNYTicks(); ++i_y) {
		double v = input_hist.GetXTick(i_y);
		std::cout << v << ", "; 
		YTicks.push_back(v);
	}
	std::cout << std::endl;
	for (int i_y = 0; i_y < input_hist.GetNYTicks(); ++i_y) std::cout << YTicks[i_y] << ", "; 
	std::cout << std::endl;

	for (int i_y = 0; i_y < input_hist.GetNYTicks()-1; ++i_y)
	for (int i_x = 0; i_x < input_hist.GetNXTicks()-1; ++i_x) {
//		Zs.push_back(input_hist.GetBinValue(i_x, i_y));
//		ZErrors.push_back(input_hist.GetBinError(i_x, i_y));
		input_hist.GetBinValue(i_x, i_y);
		input_hist.GetBinError(i_x, i_y);

		Zs.push_back(1./200.);
		ZErrors.push_back(1./200.);
	}

	std::cout << "Device iterables configuration finished. " << std::endl;
	std::cout << "XTicks size : " << XTicks.size() << std::endl;
	std::cout << "YTicks size : " << YTicks.size() << std::endl;
	std::cout << "Zs size : " << Zs.size() << std::endl;
	std::cout << "ZErrors size : " << ZErrors.size() << std::endl;

	return ArbitraryBinningHistogram2DFunctor<Iterator1, Iterator2, Iterator3, Iterator4, ArgTypeX, ArgTypeY, ArgTypeZ>(
			std::forward<Iterable1>(XTicks).begin(),
			std::forward<Iterable1>(XTicks).end(),
			std::forward<Iterable1>(YTicks).begin(),
			std::forward<Iterable1>(YTicks).end(),
			std::forward<Iterable2>(Zs).begin(),
			std::forward<Iterable2>(ZErrors).begin());
}


}  // namespace hydra


#endif /* ARBITRARYBINNINGHISTOGRAM2DFUNCTOR_H_ */
