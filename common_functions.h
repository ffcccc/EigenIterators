#ifndef __COMMON_FUNCTIONS_H__
#define __COMMON_FUNCTIONS_H__

#include <Eigen/Dense>


class CommonFunctions{
public:
	// sigmod function, depend on <cmath> library
	static double sigmod(double x) {
		return 1.0 / (1.0 + exp(-x));
	};
	
	static double crossEntropyLoss(Eigen::VectorXi y,Eigen::VectorXd h) {
		Eigen::VectorXd y_d = y.cast<double>();
		int n = y_d.size();
		double loss(0.0);
		for (int i = 0; i < n; i++) {
			loss -= (y_d(i) * log2(h(i)) + (1 - y_d(i)) * log2(1 - h(i)));
		}
		return loss / n;
	};
};



#endif