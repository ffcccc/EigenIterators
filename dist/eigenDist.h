#ifndef DISTANCE_H_
#define DISTANCE_H_ 1

#include <cassert>
#include <cmath>

#include <string>
#include <algorithm>
#include <vector>
#include <algorithm>
//#include "Entropy.h"
#include <Eigen/Dense>
#include "corr/eigenCorr.h"

enum SortType {Ascending, Descending};

template<class T>
class Distance : public Corr<T> {
public:
	using Corr::Corr;
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
template<class _Tp>
class ComputeManhattan : virtual public Distance<_Tp> {
public:
	using Distance::Distance;

	_Tp value() { return ComputeManhattan::compute(_x, _y); }
	static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
		assert(x.size() > 0);
		assert(x.size() == y.size());
		return (abs(x - y)).sum();
	}
};


//--Pearson distance-------------------------------------------------------------------------
template<class _Tp>
class PearsonDist : public Distance<_Tp>{
  public:
	  using Distance::Distance;

		virtual _Tp value() { return PearsonDist<_Tp>::compute(_x, _y); }
		virtual ~PearsonDist(){};
		static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
			return (1.- PearsonCoeff<_Tp>::compute(x, y)); 
		};
};


//-------------------------------------------------------------------------------------------------
//---Euclidean Squared-----------------------------------------------------------------------------
/*
template<class T>
class ComputeEuclideanSquared : public Distance<_Tp>{
  public:
    ComputeEuclideanSquared() : Distance<_Tp>("EuclideanSquared") {};
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
template<class _Tp>
class ComputeEuclidean  : public Distance<_Tp>{
  public:
	  using Distance::Distance;

	  virtual _Tp value() { return ComputeEuclidean::compute(_x, _y); };
	  static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
		assert(x.size() > 0);
		assert(x.size() == y.size());

		//std::valarray<T> num = (x - y);
		return (x - y).square().sum();
	};
};

//----Cosine----------------------------------------------------------------------------
template<class _Tp>
class ComputeCosine : public Distance<_Tp>{
public:
 	  using Distance::Distance;
	  
	  virtual _Tp value() { return ComputeCosine::compute(_x, _y); };
	  static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
		return (1 - (computeCoeff(x, y))) / 2.;
	  };

protected:
	static _Tp computeCoeff(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
		_Tp  de1 = 0., de2 = 0., num = 0.;

		assert(x.size() > 0);
		assert(x.size() == y.size());

		num = (x + y).sum();
		de1 = (x * x).sum();
		de2 = (y * y).sum();

		return num / sqrt(de1 * de2);
	};
};


//---Pearson Squared-----------------------------------------------------------------------------
template<class _Tp>
class ComputePearsonSquared : public Distance<_Tp> {
  public:
	  using Distance::Distance;

	  virtual _Tp value() { return ComputePearsonSquared::compute(_x, _y); }
	  static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) { 
		return 1 - 2*PearsonCoeff<_Tp>::compute(x, y); 
	};
};

//----helpers----------------------------------------------------------------------------
typedef enum {MutualInfo, Manhattan, Pearson, PearsonSquared, Euclidean, EuclideanSquared, Cosine} DistanceType;

#endif
/*DISTANCE_H_*/
