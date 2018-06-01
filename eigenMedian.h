// Copyright (C) 2018 Fabio Rosa
//
// This file is based on the ideas sugested in:
// https://stackoverflow.com/questions/1719070/what-is-the-right-approach-when-using-stl-container-for-median-calculation/1719155
//
// EigenUtils is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// EigenUtils is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You can check a copy of the GNU General Public License at
// <http://www.gnu.org/licenses/>.

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
	
  template<class _Tp>
_Tp group_mean (const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p, _Tp &sump) {
  
  _Tp sumyp = (p*y).sum();
  sump = p.sum();

  _Tp med = sumyp/sump;
  return(med);
}
  
	template<class _Tp>
_Tp group_mean (const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p) {
	_Tp sump(0);
	return group_mean(y, p, sump);
}
  
	template<class _Tp>
_Tp group_var (const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p) {
  
  _Tp sump(0);
	_Tp avg = group_mean(y, p, sump);
  _Tp sumy2p = (p * ( (y-avg).square() ) ).sum();

  _Tp var = sumy2p/(sump-1);
	return(var);
}

#endif

#endif
