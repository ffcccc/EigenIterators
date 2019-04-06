#define NOMINMAX 

#include <iostream>
#include <sstream>
#include <Eigen/Dense>
#include <limits>
#include <cstddef>

#define MyMatrix	Eigen::MatrixXd
#define MyVec		Eigen::VectorXd
#define MyType double

using Eigen::Map;
using namespace std;

#include "norm\normdens.h"
int testNorm() {
	double th(-1.334315);
	double ph(0.9089497);

	// DENSITY N(0,1)
	double t = th / ph;
	double n1 = RPORT::normdens(t, 0.0, 1.0);
	double n2 = RPORT::normdens_fast(t, 0.0, 1.0);
	double n3 = PLNK::normdens(t, 0.0, 1.0);
	double n4(0);

	n1 = RPORT::normdens(0, 0.0, 1.0);
	n2 = RPORT::normdens(0.5, 0.0, 1.0);
	n3 = RPORT::normdens(1.0, 0.0, 1.0);

	// DIST.
	// pnorm(-1.33) 										--> 0.09175914
	n1 = PLNK::pNorm01(-1.33, true);
	// pnorm(-1.33, lower.tail = FALSE) --> 0.9082409
	n2 = PLNK::pNorm01(-1.33, false);
	// pt(-1.33, 999)										--> 0.09191092
	//n3 = PLNK::pT(-1.33, 999);
	// pt(1.33, 999)										--> 0.9080891
	//n4 = PLNK::pT(1.33, 999);
	return 0;
}

#include "lr/lr.h"
bool testLR(const MyMatrix &iris){
	Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

	//data prepare,10 samples
	//Eigen::MatrixXd X(10,2);
	//X<<1.0,0.8,2.0,1.7,3.0,2.5,4.0,3.6,5.0,4.9,
	//   1.0,1.2,2.0,2.5,3.0,3.4,4.0,4.5,5.0,6.0;
	const int N(100);
	const int M(1);
	//MyVec y(N);
	MyMatrix X(N, M);
	//y = iris.col(0);
	X.col(0) = iris.col(1);

	cout << "Probs:\n" << std::defaultfloat << X.format(CleanFmt) << endl;
	Eigen::VectorXi y(10);
	//y<<0,0,0,0,0,1,1,1,1,1;
	y = iris.col(0).cast<int>();

	//train and save the weights
	LR clf1 = LR(200,0.01,0.25,0.01);  //max_iter=200,alpha=0.01(learning rate),l2_lambda=0.05,tolerance=0.01
	clf1.fit(X,y);
	cout<<"weights:\n"<<clf1.getW()<<endl; 
	//clf1.saveWeights("test.weights");

	//load the weights and predict
	LR clf2 = LR();
	//clf2.loadWeights("test.weights");
	clf2.setW(clf1.getW());
	cout << "Probs:\n" << std::defaultfloat << clf2.getW().format(CleanFmt) << endl;
	cout << "Probs:\n" <<  std::defaultfloat << clf2.predict_prob(X).format(CleanFmt) << endl;
	cout << "Predict:\n"<< std::defaultfloat << clf2.predict(X).format(CleanFmt) << endl;

	return 0;
}


#include "assoc/testAssoc.h"
bool test2x2(){
	bool res = true;
	const double PREC = 0.01;
	double OR, Conf1, Conf2, Chi2, Chi2Yates, P;
	Eigen::Matrix2d X;
	X << 2, 5, 11, 6;
	OR = computeOR(X);
	res = res && abs(OR - 0.2181818) < PREC;
	computeOR(2, 5, 11, 6, OR, Conf1, Conf2, Chi2, Chi2Yates, P);
	res = res && abs(OR - 0.2181818)<PREC;
	res = res && abs(Conf1 - 0.03205298)<PREC;
	res = res && abs(Conf2 - 1.48514424)<PREC;
	res = res && abs(P - 0.1139) < PREC;

	computeOR(34, 66, 45, 23, OR, Conf1, Conf2, Chi2, Chi2Yates, P);
	res = res && abs(OR - 0.2632997)<PREC;
	res = res && abs(Conf1 - 0.1373411)<PREC;
	res = res && abs(Conf2 - 0.5047777)<PREC;
	res = res && abs(P - 4.328e-05) < PREC;

	//computeOR(38, 82, 40, 90, OR, Conf1, Conf2, Chi2, Chi2Yates, P);
	//res = res && abs(OR - 0.87)<PREC;
	//res = res && abs(Chi2 - 0.28)<PREC;
	//res = res && abs(Chi2Yates - 0.15)<PREC;
	
	return res;
}


#include "dist/eigenDist.h"
#include "corr/eigenCorr.h"
bool testDistance(const double PREC){
	bool res = true;
	Distance<int> *distI;
	Distance<double> *distD;
	Corr<int> *corI;
	Corr<double> *corD;

	double R(0.0);

	int arrX[]	= {106,86,100,101,99,103,97,113,112,110};
	int arrY[]	= {7,0,27,50,28,29,20,12,6,17};
	double arrXd[] = { 106.,86.,100.,101.,99.,103.,97.,113.,112.,110. };
	double arrYd[] = { 7.,0.,27.,50.,28.,29.,20.,12.,6.,17. };
	
	Map<Eigen::ArrayXd> xx(&arrXd[0], 10);
	Map<Eigen::ArrayXd> yy(&arrYd[0], 10);
	Map<Eigen::ArrayXi> xi(&arrX[0], 10);
	Map<Eigen::ArrayXi>	yi(&arrY[0], 10);

	// Pearson's correlation
	 {
		R = -0.03760147;
	 	corI = new PearsonCoeff<int>(xi, yi);
	 	corD = new PearsonCoeff<double>(xx, yy);
		int Ri = corI->value();
		double Rd = corD->value();
	 	res = res && (Ri - R)	< PREC;
	 	res = res && (Rd - R)	< PREC;
		Ri = PearsonCoeff<int>::compute(xi, yi);
		Rd = PearsonCoeff<double>::compute(xx, yy);
		res = res && (Ri - R) < PREC;
		res = res && (Rd - R) < PREC;
	 	delete corI;
	 	delete corD;
	 }

	// Manhattan's correlation
	{
		distD = new ComputeManhattan<double>(xx, yy);
		double Rd = distD->value();
		res = res && (Rd - 12.0)	< PREC;
		res = res && (ComputeManhattan<double>::compute(xx, yy) - 12.0) < PREC;
		delete distD;
	}

	// Spearman's correlation
	{
		R = 0.17575757;
		ComputeSpearman<double> *distDS = new ComputeSpearman<double>(xx, yy);
		double Rd = distDS->value();
		res = res && (Rd - R) < PREC;
		res = res && (ComputeSpearman<double>::compute(xx, yy) - R) < PREC;
		delete distDS;
	}

	// Euclidean's correlation
	{
		//R =
		ComputeEuclidean<double> *distDS = new ComputeEuclidean<double>(xx, yy);
		double Rd = distDS->value();
		res = res && (Rd - R) < PREC;
		res = res && (ComputeEuclidean<double>::compute(xx, yy) - R) < PREC;
		
		delete distDS;
	}
	
	// Gamma's correlation
	{
		R = 4.44;
		ComputeGamma<double> *distDS = new ComputeGamma<double>(xx, yy);
		double Rd = distDS->value();
		res = res && (Rd - 8.333) < PREC;
		res = res && (ComputeGamma<double>::compute(xx, yy) - 8.333) < PREC;
		delete distDS;
	}
	
	// Cosine's correlation
	{
		R = 4.44;
		ComputeCosine<double> *distDS = new ComputeCosine<double>(xx, yy);
		double Rd = distDS->value();
		res = res && (Rd - 8.333) < PREC;
		res = res && (ComputeCosine<double>::compute(xx, yy) - 8.333) < PREC;
		delete distDS;
	}

	// Covariance's correlation
	{
		R = 4.44;
		ComputeCovariance<double> *distDS = new ComputeCovariance<double>(xx, yy);
		double Rd = distDS->value();
		res = res && (Rd - 8.333) < PREC;
		res = res && (ComputeEuclidean<double>::compute(xx, yy) - 8.333) < PREC;
		delete distDS;
	}

	return res;
}

#include "time\mytime.h"
#include "time\timer11.h"
int testtime()
{
	Timer tmr;
	double t = tmr.elapsed();
	std::cout << t << std::endl;

	tmr.reset();
	t = tmr.elapsed();
	std::cout << t << std::endl;
	return 0;
}

#include "fast-cpp-csv-parser-master\csv.h"
bool testIrisCSV(MyMatrix &mat) {
	io::CSVReader<5> in("../iris.csv");

	in.read_header(io::ignore_extra_column, "Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width", "Species");
	MyType val_sl;
	MyType val_sw;
	MyType val_pl;
	MyType val_pw;
	MyType y;
	string species;
	int linecount(0);
	while (in.read_row(val_sl, val_sw, val_pl, val_pw, species))
	{
		// do stuff with the data
		if (species == "setosa") {
			y = 0;
		}
		else 	if (species == "versicolor") {
			y = 1;
		}
		else 	if (species == "virginica") {
			y = 2;
		}
		mat.block<1, 5>(linecount, 0) << y, val_sl, val_sw, val_pl, val_pw;
		linecount++;
	}

	return(linecount == 150);
}

#include "lm\fastLm.h"
bool testLM(const MyMatrix &iris) {
	const int N(150);
	const int M(2);
	MyVec y(N);
	MyMatrix mat(N, M);
	y = iris.col(0);
	mat.col(0).setOnes();
	mat.col(1) = iris.col(1);
	lmsol::lmres r1 = lmsol::fastLm(mat, y, lmsol::ColPivQR_t);
	//R Call:
	//lm(formula = as.numeric(iris$Species) ~iris$Sepal.Length)

	//	Residuals :
	//	Min       1Q   Median       3Q      Max
	//	- 0.96645 - 0.40838 - 0.05035  0.33676  1.73034

	//	Coefficients :
	//	Estimate Std.Error t value Pr(> | t | )
	//	(Intercept)-2.52398    0.29878 - 8.448 2.56e-14 ***
	//	iris$Sepal.Length  0.77421    0.05063  15.292 < 2e-16 ***
	//	-- -
	//	Signif.codes : 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

	//	Residual standard error : 0.5118 on 148 degrees of freedom
	//	Multiple R - squared : 0.6124, Adjusted R - squared : 0.6098
	//	F - statistic : 233.8 on 1 and 148 DF, p - value : < 2.2e-16
	return true;
}


#include "numDeriv\numderiv.h"
int main(int argc, char **argv)
{
	if (__cplusplus == 201703L) std::cout << "C++17\n";
	else if (__cplusplus == 201402L) std::cout << "C++14\n";
	else if (__cplusplus == 201103L) std::cout << "C++11\n";
	else if (__cplusplus == 199711L) std::cout << "C++98\n";
	else std::cout << "pre-standard C++\n";

	std::cout
		<< "short: " << std::dec << std::numeric_limits<short>::min()
		<< " or " << std::hex << std::showbase
		<< std::numeric_limits<short>::min() << '\n'

		<< "int: " << std::dec << std::numeric_limits<int>::min() << std::showbase
		<< " or " << std::hex << std::numeric_limits<int>::min() << '\n' << std::dec

		<< "ptrdiff_t: " << std::numeric_limits<std::ptrdiff_t>::min() << std::showbase
		<< " or " << std::hex << std::numeric_limits<std::ptrdiff_t>::min() << '\n'

		<< "float: " << std::numeric_limits<float>::min()
		<< " or " << std::hexfloat << std::numeric_limits<float>::min() << '\n'

		<< "double: " << std::defaultfloat << std::numeric_limits<double>::min()
		<< " or " << std::hexfloat << std::numeric_limits<double>::min() << '\n';
	
	const int N(150);
	const int M(5);
	MyVec y(N);
	MyMatrix mat(N, M);

	// y "Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"
	testIrisCSV(mat);
	//
	test2x2();
	//
	testDistance(0.01);
	//
	testLR(mat.block(0,0,100,2));
	//
	testLM(mat);
	//
	testtime();
	// check numderiv working
	testDeriv();

	return 0;
}
