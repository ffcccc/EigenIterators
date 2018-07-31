#pragma once

#if !defined(_utilFile_H__INCLUDED_)
#define _utilFile_H__INCLUDED_


#include <iostream>
#include <string>
#include <valarray>

#ifdef WIN32
	//	For an MFC app, use VC_EXTRALEAN, otherwise use WIN32_LEAN_AND_MEAN. Define them before including the header files
	//	#define WIN32_LEAN_AND_MEAN
	#define VC_EXTRALEAN
	#include <windows.h>
#endif

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
using namespace boost::math;



//template<class T, class V>
//class ConfusionMatrix {
//  public:
//	inline ConfusionMatrix(const string &name_)	{ m_name = name_;	 };
//    virtual inline  ~ConfusionMatrix()				{ };
//    virtual inline const string & name()	{ return m_name;};
//    virtual inline T compute(const valarray<V> &table)	{ assert(false); return T(0);};		// array della tabella 2x2  CA CO
//																							//                 Esposti [A   B]
//																							//                 Non Espo[C   D]
//	inline V A() const { return x[0]; }
//	inline V B() const { return x[1]; }
//	inline V C() const { return x[2]; }
//	inline V D() const { return x[3]; }
//	inline V N() const { return x->sum(); }
//
//  protected:
//	string m_name;
//	valarray<V> *x;
//};
//
////--Chi2------------------------------------------------------------------------------
//template<class T, class V>
//class ComputeChi2 : virtual public ConfusionMatrix<T,V>{
//  public:
//	inline ComputeChi2() : ConfusionMatrix<T,V>("Chi Square Test") {};
//    virtual inline ~ComputeChi2(){};
//    T compute(const valarray<V> &x);
//	T compute();
//};
//
//template<class T, class V> inline
//T ComputeChi2<T,V>::compute(const valarray<V> &_x) {
//	x = &_x
//	assert(x.size() == 4);
//
//	T r		= T(0);
//	T temp	= abs((A()*D())-(C()*B()));
//	r		= (temp*temp*N())/((A+C)*(A+B)*(C+D)*(B+D));
//	return r;
//}
//
////--Chi2Yates------------------------------------------------------------------------------
//template<class T, class V>
//class ComputeChi2Yates : virtual public ConfusionMatrix<T,V>{
//  public:
//	inline ComputeChi2Yates() : ConfusionMatrix<T,V>("Chi Square Test with Yates correction") {};
//    virtual inline ~ComputeChi2Yates(){};
//    T compute(const valarray<V> &x);
//	T compute();
//};
//
//template<class T, class V> inline
//T ComputeChi2Yates<T,V>::compute(const valarray<V> &_x) {
//	x = &_x
//	assert(x.size() == 4);
//
//	T r		= T(0);
//	T temp	= abs((A()*D())-(C()*B())) - N()/(T)2;
//	r		= (temp*temp*N())/((A+C)*(A+B)*(C+D)*(B+D));
//	return r;
//}
//
////--OR------------------------------------------------------------------------------
//template<class T, class V>
//class ComputeOR : virtual public ConfusionMatrix<T,V>{
//  public:
//	inline ComputeOR() : ConfusionMatrix<T,V>("Odd Ratio") {};
//    virtual inline ~ComputeOR(){};
//    T compute(const valarray<V> &x);
//	T compute();
//};
//
//template<class T, class V> inline
//T ComputeOR<T,V>::compute(const valarray<V> &_x) {
//	x = &_x
//	assert(x.size() == 4);
//
//	T r	= (A()*D()) / (C()*B());
//	return r;
//}


//					gtypeCases[j]	gtypeControls[k]	gtypeCases[k]	gtypeControls[j]
//							A				D				B					C
static inline bool computeOR(const double &A, const double &B, const double &C, const double &D, double &OR, double &Conf1, double &Conf2, double &Chi2, double &Chi2Yates, double &P){
	// AD/BC
	if((A < 1.)||(B < 1.)||(C < 1.)||(D < 1.))
		return false;
	OR = (A*D)/(B*C);
	// Errore Standard x 1.96
	double N=A+D+B+C;
	//double ES=sqrt((OR*(1.0-OR))/N)*1.96;
	double ES=sqrt((1.0/A)+(1.0/B)+(1.0/D)+(1.0/C));
	// intervallo di confidenza
	Conf1=exp(log(OR)-1.96*ES);
	Conf2=exp(log(OR)+1.96*ES);
	// Chi Quadro
	Chi2=(pow(abs((A*D)-(C*B)),2)*N)/((A+C)*(A+B)*(C+D)*(B+D));
	//double temp=abs((A*D)-(C*B));
	//Chi2=temp*temp*N)/((A+C)*(A+B)*(C+D)*(B+D));
	chi_squared dist(1);
	//double chiP = cdf(dist, Chi2);
	// P Value
	P = 1.0 - cdf(dist, Chi2);
		//double Py=boost::math::pdf(dist, Chi2Yates);
		//double Cy=boost::math::cdf(dist, Chi2Yates);
	// Chi quadro corretto Yates
	Chi2Yates=(pow(abs((A*D)-(C*B))-(N/2),2)*N)/((A+C)*(A+B)*(C+D)*(B+D));
	//double fAB = factorial<double>(A+B);
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

static inline double computeOR(const boost::numeric::ublas::matrix<double> &ATable){
	if(0){
		double a,b,c,d;
		a=ATable(0,0);
		d=ATable(1,1);
		b=ATable(0,1);
		c=ATable(1,0);
	}
	return (ATable(0,0) * ATable(1,1)) / (ATable(0,1) * ATable(1,0));
}

static inline double computeVAR(const boost::numeric::ublas::matrix<double> &ATable){
	return 1.0/ATable(0,0) + 1.0/ATable(1,1) + 1.0/ATable(0,1) + 1.0/ATable(1,0);
}

static inline bool computeZScore(const boost::numeric::ublas::matrix<double> &ATableAff,
								 const boost::numeric::ublas::matrix<double> &ATableUnaff, double &z){
	z = fabs(	(log(computeOR(ATableAff)) - log(computeOR(ATableUnaff))) /
				sqrt( computeVAR(ATableAff)+computeVAR(ATableUnaff) )
			);
	return true;
}

#endif // !defined(_CDataFile_H__INCLUDED_)
