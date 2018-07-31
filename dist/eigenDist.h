#ifndef DISTANCE_H_
#define DISTANCE_H_ 1

#include <cassert>
#include <cmath>

#include <string>
#include <valarray>
#include <algorithm>
#include <vector>
#include <algorithm>
//#include "Entropy.h"
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::Map;

typedef enum SortType {Ascending, Descending};
typedef double TRet;

template<class T>
class Distance {
  public: 
		Distance(const std::string &name_)  { m_name = name_; };
		Distance(Map<VectorXd> x, Map<VectorXd> y, const std::string name_)  { m_name = name_; };

    virtual inline  ~Distance()      	{ /*cout << " -distance destroy- ";*/};
    virtual inline const std::string & name()	{ return m_name;};

  	virtual inline TRet compute()	{ assert(false); return (T)0;}; //= 0;
  	//virtual inline void setDims(int nx_, int ny_)	{ nX = nx_; nY = ny_; };
		
  protected:
  	std::string m_name;
  	//int nX, nY;
  	//DistanceType m_type;
};


//--MutualInfo------------------------------------------------------------------------------
//template<class T>
//class ComputeMutualInfo : virtual public Distance<T,A>{
//  public:
//	inline ComputeMutualInfo() : Distance<T,A>("MutualInfo") {};
//    virtual inline ~ComputeMutualInfo(){};
//    T compute(const std::valarray<T> &x, const std::valarray<T> &y);
//    ProbBox<T,A> pb;
//};
//
//template<class T> inline
//T ComputeMutualInfo<T,A>::compute(const std::valarray<T> &x, const std::valarray<T> &y) {
//	//assert(this->nX > 0); assert(this->nY > 0);
//	assert(x.size() > 0);
//	assert(x.size() == y.size());
//	T r = T(0);
//	pb.setDims(x,y);
//	//pb.setDims(x,this->nX,y,this->nY);
//  r = pb.mutualInformation();
//  return r;
//}

//--Manhattan------------------------------------------------------------------------------
template<class T>
class ComputeManhattan : virtual public Distance<T>{
  public:
    ComputeManhattan() : Distance<T>("Manhattan") {};
    virtual ~ComputeManhattan(){};
    virtual TRet compute(const std::valarray<T> &x, const std::valarray<T> &y);
};

template<class T> inline
TRet ComputeManhattan<T>::compute(const std::valarray<T> &x, const std::valarray<T> &y) {
    assert(x.size() > 0);
    assert(x.size() == y.size());
    return (abs(x - y)).sum();
}


//--Pearson distance-------------------------------------------------------------------------
#include "../corr/eigenCorr.h"

template<class T>
class PearsonDist : public PearsonCoeff<T>{
  public:
		PearsonDist() : PearsonCoeff<T>() {};
		virtual ~PearsonDist(){};
		virtual TRet compute(Map<VectorXd> x, Map<VectorXd> y) {
			return (1.- PearsonCoeff::compute(x, y)); 
		};
};


//-------------------------------------------------------------------------------------------------
//---Euclidean Squared-----------------------------------------------------------------------------
/*
template<class T>
class ComputeEuclideanSquared : public Distance<T>{
  public:
    ComputeEuclideanSquared() : Distance<T>("EuclideanSquared") {};
    virtual ~ComputeEuclideanSquared(){};
	  virtual TRet compute(const T *x, const T *y, const int n);

    virtual TRet compute(const std::valarray<T> &x, const std::valarray<T> &y, const std::valarray<bool> &mask);
    virtual TRet compute(const std::valarray<T> &x, const std::valarray<T> &y);
};



template<class T> inline
TRet ComputeEuclideanSquared<T>::compute(const std::valarray<T> &x, const std::valarray<T> &y, const std::valarray<bool> &mask){
	return compute(x[mask], y[mask]);
};
*/


//---Euclidean-----------------------------------------------------------------------------
template<class T>
class ComputeEuclidean  : public Distance<T>{
  public:
    ComputeEuclidean() : Distance<T>("Euclidean") {};
    virtual ~ComputeEuclidean(){};
    virtual TRet compute(const std::valarray<T> &x, const std::valarray<T> &y);
};

template<class T> inline
TRet ComputeEuclidean<T>::compute(const std::valarray<T> &x, const std::valarray<T> &y){
  assert(x.size() > 0);
  assert(x.size() == y.size());

  std::valarray<T> num = (x - y);
  //return sqrt((num * num).sum());
  return (num * num).sum();
};

//----Cosine----------------------------------------------------------------------------
template<class T>
class ComputeCosine : public Distance<T>{
  public:
    ComputeCosine() : Distance<T>("Cosine") {};
    virtual ~ComputeCosine(){};
    virtual TRet compute(const std::valarray<T> &x, const std::valarray<T> &y);

  protected:
  	virtual TRet computeCoeff(const std::valarray<T> &x, const std::valarray<T> &y);
};


template<class T> inline
TRet ComputeCosine<T>::compute(const std::valarray<T> &x, const std::valarray<T> &y) {
  return (1-(computeCoeff(x, y)))/2.;
}

template<class T> inline
TRet ComputeCosine<T>::computeCoeff(const std::valarray<T> &x, const std::valarray<T> &y) {
  T  de1 = 0., de2 = 0., num = 0.;

  assert(x.size() > 0);
  assert(x.size() == y.size());

  num = (x + y).sum();
  de1 = (x * x).sum();
  de2 = (y * y).sum();

  return num / sqrt(de1 * de2);
}


//

//---Pearson Squared-----------------------------------------------------------------------------
template<class T>
class ComputePearsonSquared : public PearsonCoeff<T>{
  public:
    ComputePearsonSquared() : PearsonCoeff<T>(/*"PearsonSquared"*/) {};
    virtual ~ComputePearsonSquared(){};

    virtual TRet compute(const std::valarray<T> &x, const std::valarray<T> &y){return 1 - 2*PearsonCoeff<T>::compute(x, y); };
};

//----helpers----------------------------------------------------------------------------
typedef enum {MutualInfo, Manhattan, Pearson, PearsonSquared, Euclidean, EuclideanSquared, Cosine} DistanceType;

#endif
/*DISTANCE_H_*/
