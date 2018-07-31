#include <stdio.h>

//#include "../math/KruskalWallisTest.h"
//
//bool testKW(){
//	bool res = true;
//	const double PREC = 0.01;
//	double KW(0.0);
//	double P(0.0);	
//	
//	int ranked[]	= {6,3,12,7,8,5,2,1,9,4,10,11};
//	int group[]		= {0,1,2,0,1,2,0,1,2,1,2,2};
//	kruskalWallisTest(group, ranked, 12, P, KW);
//	res = res && abs(KW - 5.68)	<PREC;
//	res = res && abs(P - 0.0584)<PREC;
//	
//	double data[]	= {6,3,12,7,8,5,2,1,9,4,10,11};
//	kruskalWallisTest(group, data, 12, P, KW);
//	res = res && abs(KW - 5.68)	<PREC;
//	res = res && abs(P - 0.0584)<PREC;
//	return res;
//}

#include "assoc/testAssoc.h"

bool test2x2(){
	bool res = true;
	const double PREC = 0.01;
	double OR, Conf1, Conf2, Chi2, Chi2Yates, P;
	computeOR(2, 5, 11, 6, OR, Conf1, Conf2, Chi2, Chi2Yates, P);
	res = res && abs(OR - 0.1515)<PREC;
	res = res && abs(Conf1 - 0.0223)<PREC;
	res = res && abs(Conf2 - 1.0313)<PREC;
	computeOR(34, 66, 45, 23, OR, Conf1, Conf2, Chi2, Chi2Yates, P);
	res = res && abs(OR - 2.168)<PREC;
	res = res && abs(Conf1 - 1.13)<PREC;
	res = res && abs(Conf2 - 4.16)<PREC;
	//res = res && abs( - )<PREC;
	//res = res && abs( - )<PREC;
	computeOR(38, 82, 40, 90, OR, Conf1, Conf2, Chi2, Chi2Yates, P);
	res = res && abs(OR - 0.87)<PREC;
	res = res && abs(Chi2 - 0.28)<PREC;
	res = res && abs(Chi2Yates - 0.15)<PREC;
	//res = res && abs( - )<PREC;
	//res = res && abs( - )<PREC;
	
	return res;
}


#include "dist/eigenDist.h"
#include "corr/eigenCorr.h"

bool testDistance(const double PREC){
  bool res = true;

	
	Distance<int> *distI;
	Distance<double> *distD;
	double R(0.0);
	// Pearson's correlation
	{
		int arrX[]	= {1,2,3,4};
		int arrY[]	= {1,4,9,16};
		double arrXd[]	= {1.,2.,3.,4.};
		double arrYd[]	= {1.,4.,9.,16.};
		distI = new ComputePearson<int>();
		distD = new ComputePearson<double>();
    double deb = distD->compute(arrXd, arrYd, 4);
    std::cout << deb;
		res = res && (distI->compute(arrX, arrY, 4) - 0.9844)	< PREC;
		res = res && (distD->compute(arrXd, arrYd, 4) - 0.9844)	< PREC;
		delete distI;
		delete distD;
	}

	// Manhattan's correlation
	{
		int arrX[]	= {1,6,6,7,7};
		int arrY[]	= {1,1,3,3,7};
		double arrXd[]	= {1.,6.,6.,7.,7.};
		double arrYd[]	= {1.,1.,3.,3.,7.};
		distI = new ComputeManhattan<int>();
		distD = new ComputeManhattan<double>();
		res = res && (distI->compute(arrX, arrY, 5) - 12.0)	< PREC;
		res = res && (distD->compute(arrXd, arrYd, 5) - 12.0)	< PREC;
		delete distI;
		delete distD;
	}

	// Spearman's correlation
	{
		//int arrX[]	= {106,86,100,101,99,103,97,113,112,110};
		//int arrY[]	= {7,0,27,50,28,29,20,12,6,17};
		double arrXd[]	= {106.,86.,100.,101.,99.,103.,97.,113.,112.,110.};
		double arrYd[]	= {7.,0.,27.,50.,28.,29.,20.,12.,6.,17.};
		//distI = new ComputeSpearman<double, int>();
		ComputeSpearman<double> *distDS = new ComputeSpearman<double>();
		std::valarray<double> varrX(arrXd,10);
		std::valarray<double> varrY(arrYd,10);
		R = distDS->compute(varrX, varrY);
		res = res && (R + 0.17575757)	< PREC;
		//R = distD->compute(arrXd, arrYd);
		//res = res && (R + 0.17575757)	< PREC;
		//delete distI;
		delete distDS;
//#undef USE_VALARRAY
	}

	// Euclidean's correlation
	{
		int arrX[]	= {1,2,3};
		int arrY[]	= {8,3,4};
		double arrXd[]	= {1.,2.,3.};
		double arrYd[]	= {8.,3.,4.};
		distI = new ComputeEuclidean<int>();
		distD = new ComputeEuclidean<double>();
		R = distI->compute(arrX, arrY, 3);
		res = res && (R - 7.141428)	< PREC;
		R = distD->compute(arrXd, arrYd, 3);
		res = res && (R - 7.141428)	< PREC;
		delete distI;
		delete distD;
	}
	
	// Gamma's correlation
	{
		//int arrX[]	= {1,2,3,4};
		//int arrY[]	= {1,4,9,16};
		//double arrXd[]	= {1.,2.,3.,4.};
		//double arrYd[]	= {1.,4.,9.,16.};
		//distI = new ComputeGamma<double, int>();
		//distD = new ComputeGamma<double, double>();
		//res = res && (distI->compute(arrX, arrY, 4) - 0.9844)	< PREC;
		//res = res && (distD->compute(arrXd, arrYd, 4) - 0.9844)	< PREC;
		//delete distI;
		//delete distD;
	}
	
	// Cosine's correlation
	{
		//int arrX[]	= {1,2,3,4};
		//int arrY[]	= {1,4,9,16};
		//double arrXd[]	= {1.,2.,3.,4.};
		//double arrYd[]	= {1.,4.,9.,16.};
		//distI = new ComputeCosine<double, int>();
		//distD = new ComputeCosine<double, double>();
		//res = res && (distI->compute(arrX, arrY, 4) - 0.9844)	< PREC;
		//res = res && (distD->compute(arrXd, arrYd, 4) - 0.9844)	< PREC;
		//delete distI;
		//delete distD;
	}

	// Covariance's correlation
	{
		int arrX[]	= {1,2,3,4};
		int arrY[]	= {1,4,9,16};
		double arrXd[]	= {1.,2.,3.,4.};
		double arrYd[]	= {1.,4.,9.,16.};
    distI = new ComputeCovariance<int>();
		distD = new ComputeCovariance<double>();
		res = res && (distI->compute(arrX, arrY, 4) - 8.333)	< PREC;
		res = res && (distD->compute(arrXd, arrYd, 4) - 8.333)	< PREC;
		delete distD;
	}

	return res;
}


int main(int argc, char **argv)
{
	test2x2();
	testDistance();
	
	printf("hello world\n");
	return 0;
}
