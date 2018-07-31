// Copyright Fabio Rosa 2010

#if !defined(_kruskalwallis_H__INCLUDED_)
#define _kruskalwallis_H__INCLUDED_

#ifdef _MSC_VER
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4510) // default constructor could not be generated.
#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
#endif

#include <iostream>
//#include <iomanip>
#include <map>
#include <vector>
#include <algorithm>
#include "boost/math/distributions/chi_squared.hpp"

using namespace std;

//--------------------------------------------------------------------------------
// Ascending sort function
template<class T>
struct SAscendingSort2 {
	SAscendingSort2(){};
	bool operator()(pair<T, int> rpStart, pair<T, int> rpEnd){
		return (rpStart.first < rpEnd.first);
    }
};
static inline void kruskalWallisTest(
        const int *group,
        const int *rank,
        const int nData, double &P, double &KW_stat)
{
	//
	using namespace std;
	using namespace boost::math;

	// Now we can calculate and output some stats:
	//
	double N(nData);
	map<int, double> gSize;
	map<int, double> gSum;
	for(unsigned i=0; i<N; i++){
	   gSize[group[i]]++;
	   gSum[group[i]]+=rank[i];
	}

	double sumKi(0.0);
	map<int, double>::const_iterator it = gSum.begin();
	while(it != gSum.end()){
		sumKi += pow(gSum[it->first],2)/gSize[it->first];
		//sumKi += gSize[it->first] * pow(gSum[it->first]/gSize[it->first], 2) - (3*(N+1));
		++it;
	}

	// kruskal-wallis-statistic:
	KW_stat = (12/(N*(N+1))) * sumKi  - (3*(N+1)) ;
	//cout << setw(55) << left << "KW Statistic" << "=  " << KW_stat << "\n";
	//
	//cout << setw(55) << left << "Degrees of Freedom" << "=  " << v << "\n";
	// Degrees of freedom: gSize.size() - 1
	chi_squared dist(gSize.size()-1);
	P = 1.0 - cdf(dist, KW_stat);
}

static inline void kruskalWallisTest(const int *group, const double *data, const int nData, double &P, double &KW_stat){
	//
	int *rank = new int[nData];
	vector< std::pair<double, int> > d_r(nData);

	for(int i=0;i<nData;i++){
		d_r[i].first = data[i];
		d_r[i].second = i;
	}
	sort(d_r.begin(), d_r.end(), SAscendingSort2<double>());
	for(int i=0;i<nData;i++){
		rank[ d_r[i].second ] = i+1;
	}

	kruskalWallisTest(group, rank, nData, P, KW_stat);
	delete[] rank;
}

#endif