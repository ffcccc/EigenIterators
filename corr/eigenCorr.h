#ifndef CORR_H_
#define CORR_H_ 1

#include <cassert>
#include <cmath>

#include <string>
#include <algorithm>
#include <vector>
#include <algorithm>

#include <Eigen/Dense>

template<class _Tp>
class Corr {
  public: 
	Corr(const std::string &name_)  { m_name = name_; };
	Corr(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y)  {
		//m_name = ""; 
		_x = x;
		_y = y;
	};

    virtual inline  ~Corr()      	{ /*cout << " -distance destroy- ";*/};

    //virtual inline const std::string & name()	{ return m_name; };
	
	virtual inline _Tp value() { return Corr::compute(_x, _y); };

	static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
		return _Tp(1);
	};
		
  protected:
 	Eigen::Array<_Tp, -1, 1> _x, _y;
	//std::string m_name;
	//int nX, nY;
  	//DistanceType m_type;
};



//--Pearson------------------------------------------------------------------------------
template<class _Tp>
class PearsonCoeff : public Corr<_Tp>{
public:
	using Corr::Corr;
	virtual _Tp value() { return PearsonCoeff::compute(_x, _y); };
	
	static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
		int n = x.size();
		assert(n > 0);
		assert(x.size() == y.size());

		Eigen::Array<_Tp, -1, 1> vx = (x - x.mean());
		Eigen::Array<_Tp, -1, 1> vy = (y - y.mean());

		_Tp res = (vx * vy).sum() / sqrt(vx.square().sum() * vy.square().sum());
		return res;
	};

};

//----Covariance----------------------------------------------------------------------------
template<class _Tp>
class ComputeCovariance : public Corr<_Tp>{
  public:
    using Corr::Corr;
	virtual _Tp value() { return ComputeCovariance::compute(_x, _y); };
  	
	static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
			_Tp  avx = 0., avy = 0., num = 0., n = x.size();
			assert(n > 0);
			assert(x.size() == y.size());

			avx = x.mean();
			avy = y.mean();
			num = ((x - avx) * (y - avy)).sum();

			return num / (n-1.);
	};
};


//----Gamma----------------------------------------------------------------------------
template<class _Tp>
class ComputeGamma : public Corr<_Tp>{
 public:
    using Corr::Corr;
	virtual _Tp value() { return ComputeGamma::compute(_x, _y); };

	static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
		_Tp  avx = 0., avy = 0.,
		sdevx = 0., sdevy = 0.,
		num = 0., n = x.size();

		assert(n > 0);
		assert(x.size() == y.size());

		Eigen::ArrayXd vx = (x - x.mean());
		Eigen::ArrayXd vy = (y - y.mean());

		sdevx = sqrt((vx.square()).sum());
		sdevy = sqrt((vy.square()).sum());

		num = (vx * vy).sum();
		return (num / n) / (sdevx * sdevy);
	};
};


//---Spearman-----------------------------------------------------------------------------
// todo: rivedere !! confrontare con ttest sorting
// Ascending sorting function
//template<class _Tp>
//struct SAscendingSort {
//	SAscendingSort(const Eigen::Array<_Tp, -1, 1> &VecRifPar, std::valarray<int> &VecRank){
//		VecRif = new std::valarray<_Tp>(VecRifPar);
//		for(unsigned int i=0;i<VecRank.size();i++) VecRank[i]=i; 
//	}
//	bool operator()(int rpStart, int rpEnd){
//          return (*VecRif)[rpStart] < (*VecRif)[rpEnd];
//    }
//	std::valarray<_Tp> *VecRif;
//};

// Ascending sort function
template<class _Tp>
struct SAscendingSort1 {
	SAscendingSort1(){};
	bool operator()(std::pair<_Tp, int> rpStart, std::pair<_Tp, int> rpEnd){
		return (rpStart.first < rpEnd.first);
  }
};

template<class _Tp>
class ComputeSpearman : public Corr<_Tp>{
  public:
    using Corr::Corr;
	virtual _Tp value() { return ComputeSpearman::compute(_x, _y); };

    static _Tp compute(const Eigen::Array<_Tp, -1, 1> &x, const Eigen::Array<_Tp, -1, 1> &y) {
			unsigned int n = x.rows();
			assert(n > 0);
			assert(x.rows() == y.rows());
			_Tp d;
			Eigen::Array<int, -1, 1> xRank(n), yRank(n);
			std::vector< std::pair<double, int> > d_x(n);
			std::vector< std::pair<double, int> > d_y(n);
			for(int i=0;i<n;i++){
				d_x[i].first = x[i];
				d_x[i].second = i;
				d_y[i].first = y[i];
				d_y[i].second = i;
			}
			std::sort(d_x.begin(), d_x.end(), SAscendingSort1<_Tp>());
			std::sort(d_y.begin(), d_y.end(), SAscendingSort1<_Tp>());
			for(int i=0;i<n;i++){
				xRank( d_x[i].second ) = i+1;
				yRank( d_y[i].second ) = i+1;
			}

			//std::sort(&xRank[0], &xRank[0]+xRank.size(),SAscendingSort<_Tp>(x,xRank));
			//std::sort(&yRank[0], &yRank[0]+yRank.size(),SAscendingSort<_Tp>(y,yRank));
				
			//std::valarray<int> num = (xRank - yRank);
			d = (xRank - yRank).square().sum();
			//return (1-(1 - (6*d / (n*(n*n-1)))))/2;
			
			// todo: check
			return 1.0 - 6.0 * d / (n*(n*n-1));
			//return (3*d / (n*(n*n-1)));
	};
};
//


//----helpers----------------------------------------------------------------------------
typedef enum {PearsonC, Covariance, Spearman, Gamma/*, MCI, HGG, dcor*/} CorrType;

#endif
/*CORR_H_*/
