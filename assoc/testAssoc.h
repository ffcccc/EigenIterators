#pragma once

#if !defined(_utilFile_H__INCLUDED_)
#define _utilFile_H__INCLUDED_


#include <iostream>
#include <string>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace std;
using namespace boost::math;


//					gtypeCases[j]	gtypeControls[k]	gtypeCases[k]	gtypeControls[j]
//							A				D				B					C
static inline bool computeOR(const double &A, const double &B, const double &C, const double &D, double &OR, double &Conf1, double &Conf2, double &Chi2, double &Chi2Yates, double &P){
	// AD/BC
	if((A < 1.)||(B < 1.)||(C < 1.)||(D < 1.))
		return false;

	// OR
	OR = (A*D)/(B*C);

	// Errore Standard
	double N=A+D+B+C;
	//double ES=sqrt((OR*(1.0-OR))/N)*1.96;
	double ES=sqrt((1.0/A)+(1.0/B)+(1.0/D)+(1.0/C));

	// intervallo di confidenza: Errore Standard x 1.96
	Conf1=exp(log(OR)-1.96*ES);
	Conf2=exp(log(OR)+1.96*ES);

	// Chi Quadro
	Chi2=(pow(abs((A*D)-(C*B)),2)*N)/((A+C)*(A+B)*(C+D)*(B+D));
	//double temp=abs((A*D)-(C*B));
	//Chi2=temp*temp*N)/((A+C)*(A+B)*(C+D)*(B+D));

	// P Value
	chi_squared dist(1);
	P = 1.0 - cdf(dist, Chi2);

	// Chi quadro corretto Yates
	Chi2Yates=(pow(abs((A*D)-(C*B))-(N/2),2)*N)/((A+C)*(A+B)*(C+D)*(B+D));

	// P di Fischer
	//double pFischer=(factorial<double>(A+B)*factorial<double>(C+D)*factorial<double>(A+C)*factorial<double>(B+D))
	//	 /(factorial<double>(N)*factorial<double>(A)*factorial<double>(B)*factorial<double>(C)*factorial<double>(D));
	return true;
}

//
static inline bool computeOR(const double &A, const double &B, const double &C, const double &D, double &OR, double &Chi2, double &P){
	// AD/BC
	bool doYates=false;
	bool lowFreq=false;
	Chi2=1.;
	if((A < 5)||(B < 5)||(C < 5)||(D < 5)) 
		lowFreq=true;
	if(lowFreq && ((A < 1.)||(B < 1.)||(C < 1.)||(D < 1.))) 
		return false;
	OR = (A*D)/(B*C);
	
	double N=A+D+B+C;
	if((N < 20) && lowFreq){
		// Chi quadro corretto Yates	
		doYates=true;
		Chi2=(pow(abs((A*D)-(C*B))-(N/2),2)*N)/((A+C)*(A+B)*(C+D)*(B+D));
	} else {
		// Chi Quadro
		Chi2=(pow(abs((A*D)-(C*B)),2)*N)/((A+C)*(A+B)*(C+D)*(B+D));
	}
	// P Value
	chi_squared dist(1);
	P = 1.0 - cdf(dist, Chi2);
	return true;
}

static inline bool computeClassErr(const double &A, const double &B, const double &C, const double &D, double &BA){
	if((A < 1.)||(B < 1.)||(C < 1.)||(D < 1.))
		return false;
	BA = (A/(A+C) + D/(D+B))/2.0;
	return true;
}

template<class _Tp>
static inline double computeOR(const Eigen::Matrix<_Tp, 2, 2> &ATable){
	//if(0){
	//	double a,b,c,d;
	//	a=ATable(0,0);
	//	d=ATable(1,1);
	//	b=ATable(0,1);
	//	c=ATable(1,0);
	//}
	return (ATable(0,0) * ATable(1,1)) / (ATable(0,1) * ATable(1,0));
}

template<class _Tp>
static inline double computeVAR(const Eigen::Matrix<_Tp, 2, 2> &ATable){
	_Tp result;
	//result = 1.0 / ATable(0, 0) + 1.0 / ATable(1, 1) + 1.0 / ATable(0, 1) + 1.0 / ATable(1, 0);
	result = ATable.cwiseInverse().sum();
	return result;
}

template<class _Tp>
static inline bool computeZScore(const Eigen::Matrix<_Tp, 2, 2> &ATableAff,
								 const Eigen::Matrix<_Tp, 2, 2> &ATableUnaff, double &z){
	z = fabs(	(log(computeOR<_Tp>(ATableAff)) - log(computeOR<_Tp>(ATableUnaff))) /
				sqrt( computeVAR<_Tp>(ATableAff)+computeVAR<_Tp>(ATableUnaff) )
			);
	return true;
}

#endif // !defined(_CDataFile_H__INCLUDED_)
