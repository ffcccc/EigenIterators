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

//#include <random>
#include <algorithm>
#include <vector>
#include "iter/eigenBeginEnd.h"

// MSVC tip: use /Zc:__cplusplus compiler language option
#if __cplusplus >= 201103L


template<class _Tp>
inline _Tp var(const Eigen::Array<_Tp, -1, 1>& v)
{
	_Tp N(v.size());
	if (N < 2) {
		return _Tp(0);
	}

	_Tp v_avg = v.mean();
	_Tp res = (v - v_avg).square().sum();
	return res / (N - 1);
}

template<class _Tp>
inline _Tp stdDev(const Eigen::Array<_Tp, -1, 1>& v)
{
	_Tp res = std::sqrt(var(v));
	return res;
}

// todo: test it !!!
// to verify implementation: https://stackoverflow.com/questions/15138634/eigen-is-there-an-inbuilt-way-to-calculate-sample-covariance#15142446
template<class _Tp>
inline _Tp cov(const Eigen::Matrix<_Tp, -1, -1>& mat, Eigen::Matrix<_Tp, -1, -1>& matcov)
{
	_Tp N(mat.nrows());
	_Tp M(mat.ncols());
	if ((N < 2) || (M < 2)) {
		return _Tp(0);
	}

	Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
	matcov = (centered.adjoint() * centered) / double(mat.rows() - 1);

	return _Tp(1);
}

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
	if (v.empty()) {
		return _Tp(0);
	}
	auto n = v.size() / 2;
	std::nth_element(v.begin(), v.begin() + n, v.end());
	_Tp med = v[n];
	if (!(v.size() & 1)) { //If the set size is even
		_Tp max_it = *std::max_element(v.begin(), v.begin() + n);
		med = (max_it + med) / 2.0;
	}
	return med;
}


template<class _Tp>
inline _Tp median(Eigen::Array<_Tp, -1, 1>& v)
{
	switch (v.size()) {
	case 0:
		// assert ?
		return _Tp(0);
		break;
	case 1:
		return v(0);
		break;
	case 2:
		return v.mean();
		break;
	}
	auto n = v.size() / 2;
	std::nth_element(begin(v), begin(v) + n, end(v));
	_Tp med = v(n);
	if (!(v.size() & 1)) { //If the set size is even
		_Tp max_it = v.head(n).maxCoeff();
		med = (max_it + med) / 2.0;
	}
	return med;
}

// todo:
// it works only for positive val of y
template<class _Tp>
inline _Tp group_median(Eigen::Array<_Tp, -1, 1>& v, const Eigen::Array<_Tp, -1, 1>&  p)
{
	if (v.size() == 0) {
		return _Tp(0);
	}
	_Tp yp = (p*v);
	int n_yp = p.sum();
	int n_exclude = v.size() - n_yp;

	auto n = n_yp / 2 + n_exclude;
	std::nth_element(begin(v), begin(v) + n, end(v));
	_Tp med = v(n);
	if (!(v.size() & 1)) { //If the set size is even
		_Tp max_it = v.head(n).maxCoeff();
		med = (max_it + med) / 2.0;
	}
	return med;
}

template<class _Tp>
inline _Tp group_mean(const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p, _Tp &sump) {

	_Tp sumyp = (p*y).sum();
	sump = p.sum();

	_Tp med = sumyp / sump;
	return(med);
}

template<class _Tp>
inline _Tp group_mean(const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p) {
	_Tp sump(0);
	return group_mean(y, p, sump);
}

template<class _Tp>
inline _Tp group_var(const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p) {

	_Tp sump(0);
	_Tp avg = group_mean(y, p, sump);
	_Tp sumy2p = (p * ((y - avg).square())).sum();

	_Tp var = sumy2p / (sump - 1);
	return(var);
}

template<class _Tp>
inline _Tp group_stdDev(const Eigen::Array<_Tp, -1, 1>& y, const Eigen::Array<_Tp, -1, 1>&  p) {
	_Tp sump(0);
	return std::sqrt(group_var(y, p, sump));
}

//
// Sm1 = Sample Mean 1.
// Sd1 = Sample Standard Deviation 1.
// Sn1 = Sample Size 1.
// Sm2 = Sample Mean 2.
// Sd2 = Sample Standard Deviation 2.
// Sn2 = Sample Size 2.
// alpha = Significance Level.
//
// A Students t test applied to two sets of data.
// We are testing the null hypothesis that the two
// samples have the same mean and that any difference
// if due to chance.
// See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
#include <boost/math/distributions/students_t.hpp>
template<class _Tp>
inline _Tp two_samples_t_test_equal_sd(
	_Tp Sm1,
	_Tp Sd1,
	unsigned Sn1,
	_Tp Sm2,
	_Tp Sd2,
	unsigned Sn2,
	_Tp alpha)
{

	//
	using namespace std;
	using namespace boost::math;

	// Now we can calculate and output some stats:
	//
	// Degrees of freedom:
	_Tp v = Sn1 + Sn2 - 2;
	//cout << setw(55) << left << "Degrees of Freedom" << "=  " << v << "\n";
	// Pooled variance:
	_Tp sp = std::sqrt(((Sn1 - 1) * Sd1 * Sd1 + (Sn2 - 1) * Sd2 * Sd2) / v);
	//cout << setw(55) << left << "Pooled Standard Deviation" << "=  " << v << "\n";
	// t-statistic:
	_Tp t_stat = (Sm1 - Sm2) / (sp * sqrt(1.0 / Sn1 + 1.0 / Sn2));
	students_t dist(v);
	_Tp q = cdf(complement(dist, fabs(t_stat)));
	return q;
}

template<class _Tp>
inline _Tp two_samples_t_test_unequal_sd(
	double Sm1,
	double Sd1,
	unsigned Sn1,
	double Sm2,
	double Sd2,
	unsigned Sn2,
	double alpha)
{
	//
	// Sm1 = Sample Mean 1.
	// Sd1 = Sample Standard Deviation 1.
	// Sn1 = Sample Size 1.
	// Sm2 = Sample Mean 2.
	// Sd2 = Sample Standard Deviation 2.
	// Sn2 = Sample Size 2.
	// alpha = Significance Level.
	//
	// A Students t test applied to two sets of data.
	// We are testing the null hypothesis that the two
	// samples have the same mean and that any difference
	// if due to chance.
	// See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
	//
	using namespace std;
	using namespace boost::math;

	// Now we can calculate and output some stats:
	//
	// Degrees of freedom:
	_Tp v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;
	v *= v;
	_Tp t1 = Sd1 * Sd1 / Sn1;
	t1 *= t1;
	t1 /= (Sn1 - 1);
	_Tp t2 = Sd2 * Sd2 / Sn2;
	t2 *= t2;
	t2 /= (Sn2 - 1);
	v /= (t1 + t2);
	// t-statistic:
	_Tp t_stat = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);

	// Define our distribution, and get the probability:
	students_t dist(v);
	_Tp q = cdf(complement(dist, fabs(t_stat)));
	return q;
}

#endif

//#include <bitset>
#include <list>
// sample(x, y, replace = FALSE, prob = NULL)
//
//	x				a vector of one or more elements from which to choose
//	y				a vector of the desired size of elements uniformly choosen from x
//	replace	should sampling be with replacement?
//	not yet impl:
//	prob		a vector of probability weights for obtaining the elements of the vector being sampled.
template<class _Tp>
inline bool sample(const Eigen::Array<_Tp, -1, 1>& x, Eigen::Array<_Tp, -1, 1>& y, bool replace = false, bool use_seed = true, unsigned int seed = 12345)
{
	int sn = x.size();	// uniform max param
	int sk = y.size(); // num items to choose
	if (sk == 0 || sn == 0) {
		return false;
	}

	// construct a marsenne twister random generator engine from a user-seed or randomly choosen OS request (chrono or whatever...):
	unsigned int localseed(seed);
	if (!use_seed) {
		mutable std::random_device rd;
		localseed = rd();
		//localseed = 18;
	}

	mutable std::mt19937 mt(localseed);

	mutable std::uniform_int_distribution dist(0, sn - 1);
	if (replace) {
		for (int i = 0; i < sk; ++i) {
			auto x_pos = dist(mt);
			y(i) = x(x_pos);
		}
	}
	// no replace
	else {
		std::list<double> xl(begin(x), end(x));
		std::list<double>::iterator itdel;
		for (int i = 0; i < sk; ++i) {
			itdel = xl.begin();
			mutable std::uniform_int_distribution::param_type parms(0, xl.size() - 1);
			dist.param(parms);
			auto x_pos = dist(mt);
			y(i) = x(x_pos);
			// now remove already choosen items
			advance(itdel, x_pos);
			xl.erase(itdel);

		}
	}

	return true;
}

// sample(x, size, replace = FALSE, prob = NULL)
//
//	x				either a vector of one or more elements from which to choose, or a positive integer. See ‘Details.’
//	size		a non-negative integer giving the number of items to choose.
//	replace	should sampling be with replacement?
//	not yet impl:
//	prob		a vector of probability weights for obtaining the elements of the vector being sampled.
template<class _Tp>
inline Eigen::Array<_Tp, -1, 1> sample(const Eigen::Array<_Tp, -1, 1>& x, const unsigned int ysize, bool replace = false, bool use_seed = true, unsigned int seed = 12345)
{
	Eigen::Array<_Tp, -1, 1> y(ysize);
	bool res = sample(x, y, replace, use_seed, seed);
	if (res)
		return y;

	y.setConstant(ysize, _Tp(0));
	return y;
}


#endif
