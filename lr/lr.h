#ifndef __LR_H__
#define __LR_H__

#include "common_functions.h"
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class LR{
public:
	LR(int max_i=100,double alp=0.01,double l2_lambda=0.05,double tolerance=0.01) //the params name can't be the same as the class member?
	{
		lambda = l2_lambda;
		max_iter = max_i;
		tol = tolerance;
		alpha = alp;
	}

	~LR() {};
	void fit(Eigen::MatrixXd X,Eigen::VectorXi y) {
		//learn VectorXd W, consider reg,max_iter,tol.   
		//TODO: check X,y

		//VectorXd W = Eigen::VectorXd::Random(X.cols()+1); wrong! u can not declare W again,otherwise it didn't represent the class member
		W = Eigen::VectorXd::Random(X.cols() + 1);  //the last column of weight represent b
		Eigen::MatrixXd X_new(X.rows(), X.cols() + 1);
		X_new << X, Eigen::MatrixXd::Ones(X.rows(), 1);  //last column is 1.0
		double loss;
		Eigen::VectorXd E;
		for (int iter = 0; iter < max_iter; iter++) {
			Eigen::VectorXd y_pred = predict_prob(X);
			Eigen::VectorXd y_d = y.cast<double>();  //cast type first
			E = y_pred - y_d;

			//W = (1.0-lambda/y.size())*W - alpha*X_new.transpose()*E;  //W:= (1-lambda/n_samples)W-alpha*X^T*E
			Eigen::VectorXd dw = (E.transpose()*X_new) / y.size();
			W = W - lambda * dw;
			//reference : http://blog.csdn.net/pakko/article/details/37878837

			//when loss<tol, break
			loss = CommonFunctions::crossEntropyLoss(y, predict_prob(X));
			if (loss <= tol)
				break;
		}

	};

	Eigen::VectorXd getW() {
		return W;
	};

	Eigen::Index setW(const Eigen::VectorXd &_w) {
		W = _w;
		return W.size();
	};

	Eigen::VectorXd predict_prob(Eigen::MatrixXd X) {
		//predict the probability (of label 1) for given data X
		Eigen::MatrixXd X_new(X.rows(), X.cols() + 1);
		X_new << X, Eigen::MatrixXd::Ones(X.rows(), 1);
		int num_samples = X_new.rows();
		Eigen::VectorXd y_pred_prob = Eigen::VectorXd::Zero(num_samples);
		for (int num = 0; num < num_samples; num++) {
			y_pred_prob(num) = CommonFunctions::sigmod(X_new.row(num).dot(W));
		}

		return y_pred_prob;
	};
	Eigen::VectorXi predict(Eigen::MatrixXd X) {
		//predict the label for given data X
		Eigen::VectorXd y_pred_prob = predict_prob(X);
		Eigen::VectorXi y_pred(y_pred_prob.size());
		for (int num = 0; num < y_pred_prob.size(); num++) {
			y_pred(num) = y_pred_prob(num) > 0.5 ? 1 : 0;
		}
		return y_pred;
	};
	void saveWeights(std::string filename) {
		//save the model (save the weight ) into filename. 
		std::ofstream ofile;
		std::string path = "./weights/" + filename;
		ofile.open(path.c_str());
		if (!ofile.is_open()) {
			std::cerr << "Can not open the file when call LR::saveParams" << std::endl;
			return;
		}

		//W write into the file
		for (int i = 0; i < W.size() - 1; i++) {
			ofile << W(i) << " ";
		}
		ofile << W(W.size() - 1);
		ofile.close();
	};
	void loadWeights(std::string filename) {
		//load the model (load the weight ) from filename.
		std::ifstream ifile;
		std::string path = "./weights/" + filename;
		ifile.open(path.c_str());
		if (!ifile.is_open()) {
			std::cerr << "Can not open the file when call LR::loadParams" << std::endl;
			return;
		}

		//read the weights into vector<double>
		std::string line;
		std::vector<double> weights;
		getline(ifile, line);    //only one line
		std::stringstream ss(line);
		double tmp;
		while (!ss.eof()) {
			ss >> tmp;
			weights.push_back(tmp);
		}

		//initialize VectorXd with std::vector
		W = Eigen::VectorXd::Map(weights.data(), weights.size());

		ifile.close();
	};
private:
	Eigen::VectorXd W;
	int max_iter;
	double lambda;  //l2 regulization
	double tol;  // error tolence
	double alpha;
};



#endif