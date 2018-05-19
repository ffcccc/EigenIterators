#ifndef ___EigenHelper___
#define ___EigenHelper___

#include <algorithm>
#include <vector>
#include "eigenBeginEnd.h"

#if __cplusplus >= 201103L
  /**
   *  @brief  post-condition: After returning, the elements in v may be reordered and the resulting order is implementation defined.
   * 		https://stackoverflow.com/questions/1719070/what-is-the-right-approach-when-using-stl-container-for-median-calculation/1719155#1719155         
   *		This algorithm handles both even and odd sized inputs efficiently using the STL nth_element (amortized O(N)) algorithm
	*		and the max_element algorithm (O(n)). Note that nth_element has another guaranteed side effect, namely that all of the elements
	*		before n are all guaranteed to be less than v[n], just not necessarily sorted.
   *  @param  __va  valarray.
   */
  template<class _Tp>
    inline _Tp median(std::vector<_Tp> &v)
	{
	  if(v.empty()) {
		return _Tp(0);
	  }
	  auto n = v.size() / 2;
	  std::nth_element(v.begin(), v.begin()+n, v.end());
	  _Tp med = v[n];
	  if(!(v.size() & 1)) { //If the set size is even
		_Tp max_it = *std::max_element(v.begin(), v.begin()+n);
		med = (max_it + med) / 2.0;
	  }
	  return med;    
	}


  template<class _Tp>
    inline _Tp median(Eigen::Array<_Tp, -1, 1>& v)
	{
	  if(v.size() == 0) {
		return _Tp(0);
	  }
	  auto n = v.size() / 2;
	  std::nth_element(begin(v), begin(v)+n, end(v));
	  _Tp med = v(n);
	  if(!(v.size() & 1)) { //If the set size is even
		_Tp max_it = v.head(n).maxCoeff();
		med = (max_it + med) / 2.0;
	  }
	  return med;    
	}
	
//	double fuzzy_median (std::vector<double> y, std::vector<double> p) {
//		std::vector<MyType> groups[3];
//		for(auto i=0; i<N; i++){
//			int pos = int(snp(i));
//			groups[pos].push_back(y(i));
//		}
//		for(auto i=0; i<3; i++){
//			med[i] = median(groups[i]);
//		}
//
//		return(med);
//	}
	
	double fuzzy_mean (std::vector<double> y, std::vector<double> p) {
  
  int n = p.size();

  double sumyp = 0;
  double sump = 0;
  for (int i=0;i<n;++i){ 
    sumyp += p[i]*y[i];
    sump += p[i];
	}

  double med = sumyp/sump;

  return(med);
}

  template<class _Tp>
_Tp fuzzy_mean (const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p, _Tp &sump) {
  
  _Tp sumyp = (p*y).sum();
  sump = p.sum();

  _Tp med = sumyp/sump;
  return(med);
}
  
	template<class _Tp>
_Tp fuzzy_mean (const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p) {
	_Tp sump(0);
	return fuzzy_mean(y, p, sump);
}
  
	template<class _Tp>
_Tp fuzzy_var (const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p) {
  
  _Tp sump(0);
	_Tp avg = fuzzy_mean(y, p, sump);
  _Tp sumy2p = (p * ( (y-avg).square() ) ).sum();

  _Tp var = sumy2p/(sump-1);
	return(var);
}

#endif

#endif
