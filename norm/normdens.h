#pragma once

#ifndef ___NORMDENS___
#define ___NORMDENS___

//
// this module was developed with version 3.3.4 of Eigen Matrix lib
// it refers the Eigen lib path through the link Eigen334 to avoid 
// issues when building the corresponding R package Remcpp
//

#include "Eigen/Core"
#include "Eigen/Dense"

// header EM
//#include <math.h>
#ifdef EM_USEFLOAT
#define MyMatrix	Eigen::MatrixXf
#define MyVec		Eigen::VectorXf
#define MyType float
#else
#define MyMatrix	Eigen::MatrixXd
#define MyVec		Eigen::VectorXd
#define MyType double
#endif

namespace RPORT {
	// same as R pi: > sprintf("%.100f",pi)
	const MyType _PI(3.1415932653589793115997963468544185161590576171875);
	// R/.../src/include/Rmath.h0.in
	// from abramowitz & stegun, p. 3.
	const MyType M_1_SQRT_2PI(0.398942280401432677939946059934);	// 1/sqrt(2pi)
	//	const MyType M_LOG10E(0.434294481903251827651128918917);	/* log10(e) */
	const MyType M_LN2(0.693147180559945309417232121458);	/* ln(2) */
//	const MyType M_LN10(2.302585092994045684017991454684);	/* ln(10) */

	//MyVec normdens(MyVec x, MyVec muu, MyType sig = 1.0);
	//MyType normdens(MyType x, MyType muu = 0.0, MyType sig = 1.0);
	//
	//MyVec normdens_fast(MyVec x, MyVec muu, MyType sig = 1.0);
	//MyType normdens_fast(MyType x, MyType muu = 0.0, MyType sig = 1.0);

//--------------------------------------------------------------------------
// Normal PDF: R technique
//
// from  r-source/src/nmath/dnorm.c MATHLIB_FAST_dnorm
	inline MyType normdens(MyType v, MyType mu, MyType sigma)
	{
		MyType res(0.);
		MyType x = (v - mu) / sigma;
		x = fabs(x);
		if (x < 5) {
			res = RPORT::M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
		}
		else {
			if (x > sqrt(-2 * M_LN2*(DBL_MIN_EXP + 1 - DBL_MANT_DIG))) {
				res = 0.;
			}
			else {
				double x1 = ldexp(round(ldexp(x, 16)), -16);
				double x2 = x - x1;
				res = M_1_SQRT_2PI / sigma * (exp(-0.5 * x1 * x1) * exp((-0.5 * x2 - x1) * x2));
			}
		}
		return res;
	};

	inline MyVec normdens(MyVec v, MyVec mu, MyType sigma)
	{
		MyVec res = v;
		for (auto i = 0; i < v.size(); i++) {
			res(i) = normdens(v(i), mu(i), sigma);
		}
		return res;
	};

	//--------------------------------------------------------------------------
	// Normal PDF: R technique
	//
	// from  r-source/src/nmath/dnorm.c MATHLIB_FAST_dnorm
	inline MyType normdens_fast(MyType v, MyType mu, MyType sigma)
	{
		MyType x = (v - mu) / sigma;
		x = fabs(x);
		MyType res = RPORT::M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
		return res;
	};

	inline MyVec normdens_fast(MyVec v, MyVec mu, MyType sigma)
	{
		MyVec x = -(v.array() - mu.array()) / sigma;
		// x = fabs(x)  --> skipped due to ^2
		MyVec res = (RPORT::M_1_SQRT_2PI / sigma) * exp(-0.5 * x.array().square());
		return res;
	};
}

//--------------------------------------------------------------------------
// Normal PDF: plink approach

namespace PLNK {
	const MyType _PI(3.141593);
	//MyVec normdens(MyVec x, MyVec muu, MyType sig = 1.0);
	//MyVec normdens(MyVec x, MyType muu = 0.0, MyType sig = 1.0);
	//MyType normdens(MyType x, MyType muu = 0.0, MyType sig = 1.0);

	inline MyType normdens(MyType x, MyType muu, MyType sig)
	{
		MyType sig2 = 2 * sig * sig;
		MyType res = 1 / (sqrt(PLNK::_PI * sig2)) * exp(-pow(x - muu, 2) / sig2);
		return res;
	};

	// x:vector, muu:array of means, sigma: const std
	inline MyVec normdens(MyVec x, MyVec muu, MyType sig)
	{
		MyType sig2 = 2 * sig * sig;
		MyType fact1 = 1 / (sqrt(PLNK::_PI * sig2));

		MyVec fact2 = -(x.array() - muu.array()).square() / sig2; //* coef1;
		MyVec res = fact1 * exp(fact2.array());
		return res;
	};

	// x:vector, muu:const mean, sigma: const std
	inline MyVec normdens(MyVec x, MyType muu, MyType sig)
	{
		MyType sig2 = 2 * sig * sig;
		MyType fact1 = 1 / (sqrt(PLNK::_PI * sig2));

		MyVec fact2 = -(x.array() - muu).square() / sig2; //* coef1;
		MyVec res = fact1 * exp(fact2.array());
		return res;
	};


	inline double pNorm01(double z, bool lowertail)
	{
		double sqrt2pi = 2.50662827463;
		double t0, z1, p0;
		t0 = 1 / (1 + 0.2316419 * fabs(z));
		z1 = exp(-0.5 * z*z) / sqrt2pi;
		p0 = z1 * t0
			* (0.31938153 +
				t0 * (-0.356563782 +
					t0 * (1.781477937 +
						t0 * (-1.821255978 +
							1.330274429 * t0))));
		if (lowertail) {
			return z >= 0 ? 1 - p0 : p0;
		}
		return z >= 0 ? p0 : 1 - p0;
	}

}

#endif