// Copyright (C) 2007 Fabio Rosa rosa@isi.it
//
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//    derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Entropy.h"
#include <iostream>

// http://www.onlamp.com/pub/a/php/2005/03/24/joint_entropy.html?page=1

//	 Y	buys_computer 	 
// X	no	yes	|	Σi+
//age0 	3 	2 	|	5
//age1 	0 	4 	|	4
//age2 	2 	3 	|	5
//--------------------
//Σ+j	5 	9 	|	14

//	 Y	buys_computer 	 
// X	no 			yes 	|	Σi+
//age0 	0.21429 	0.14286 |	0.35714
//age1	0.00000 	0.28571 |	0.28571
//age2 	0.14286 	0.21429 |	0.35714
//--------------------------------------
//Σ+j 	0.35714 	0.64286 |	1

//H(X,Y) = -1 * [ (0.21429 * log(0.21429)) + 
//                (0.14286 * log(0.14286)) + 
//                (0.28571 * log(0.28571)) + 
//                (0.14286 * log(0.14286)) + 
//                (0.21429 * log(0.21429)) ]
//
//H(X,Y) = -1 * [ -0.476226947429 + 
//                -0.401050703151 + 
//                -0.516387120588 + 
//                -0.476226947429 + 
//                -0.401050703151 ]
//
//H(X,Y) = -1 * -2.270942421748
//
//H(X,Y) = 2.271
//H(buys_computer) = 0.940285958671
//H(buys_computer | age) = 0.693536138896
//I(age;buys_computer) = H(buys_computer) - H(buys_computer | age) = 0.246749819774

float abuf[13] = {0,0,1,1,1,1,0,0,0,1,0,1,1};
float bbuf[13] = {0,0,1,1,1,0,1,0,1,1,1,1,1};

void entropyUnitTest(double eps){
	float agebuf[14] = {0,0,1,2,2,2,1,0,0,2,0,1,1,2};
	valarray<float> age(agebuf, 14);
	float buybuf[14] = {0,0,1,1,1,0,1,0,1,1,1,1,1,0};
	valarray<float> buy(buybuf, 14);	
//	valarray<int> va1 ( 6 );
//		va1 [ 0 ] = 3;  va1 [ 1 ] = 0;  va1 [ 2 ] = 2;
//		va1 [ 3 ] = 2;  va1 [ 4 ] = 4;  va1 [ 5 ] = 3;
//  double py0=5,  py1=9,  px0=5,  px1=4,  px2=5;
  
  ProbBox<float,float> *pb = new ProbBox<float,float>(age,buy);	//age,3,buy,2);
  
  cout << "P(Y=1) 0.64286 -> "<< pb->prob(Y,1) << std::endl;// assert(< eps);
  cout << "P(Y=0) 0.35714 -> "<< pb->prob(Y,0) << std::endl;// assert(< eps); 
  cout << "P(X=0) 0.35714 -> "<< pb->prob(X,0) << std::endl;// assert(< eps);
  cout << "P(X=0) 0.28571 -> "<< pb->prob(X,1) << std::endl;// assert(< eps);
  cout << "P(X=0) 0.35714 -> "<< pb->prob(X,2) << std::endl;// assert(< eps);
  cout << std::endl;
  cout << "P(Y=0|X=0) 0.6 -> "<< pb->conditionalProbYX(0,0) << std::endl;// assert(< eps); 
  cout << "P(Y=0|X=1) 0.0 -> "<< pb->conditionalProbYX(1,0) << std::endl;// assert(< eps); 
  cout << "P(Y=0|X=2) 0.4 -> "<< pb->conditionalProbYX(2,0) << std::endl;// assert(< eps);
  cout << "P(Y=1|X=0) 0.4 -> "<< pb->conditionalProbYX(0,1) << std::endl;// assert(< eps); 
  cout << "P(Y=1|X=1) 1.0 -> "<< pb->conditionalProbYX(1,1) << std::endl;// assert(< eps); 
  cout << "P(Y=1|X=2) 0.6 -> "<< pb->conditionalProbYX(2,1) << std::endl;// assert(< eps);
  cout << std::endl;
  cout << "H(X,Y) = 2.271 -> " << pb->jointEntropyXY() << std::endl;// assert(< eps);
  cout << "H(buys_computer) = 0.940285958671 -> " << pb->entropy(Y) << std::endl;// assert(< eps);
  cout << "H(age) =  -> " << pb->entropy(X) << std::endl;// assert(< eps);
  cout << "I(age;buys_computer) = H(buys_computer) - H(buys_computer | age) = 0.246749819774 -> " << pb->mutualInformation()  << std::endl;// assert(< eps);
  delete pb;
  
  ProbBox<float,float> pb2;
//  for(int i=0; i<3; i++){
  valarray<float> a(abuf, 13);
  valarray<float> b(bbuf, 13);	
  	
  //pb2.setDims(age,buy);	//age,3,buy,2);
  pb2.updateContingencyTable(age,buy);
  float r = pb2.mutualInformation();
  cout << "r: " << r << std::endl;
  //pb2.setDims(buy,age);
  pb2.updateContingencyTable(buy,age);
  float r2 = pb2.mutualInformation();
  cout << "r: " << r2 << std::endl;
  //pb2.setDims(b,a);
  pb2.updateContingencyTable(a,b);
  float r3 = pb2.mutualInformation();
  cout << "r: " << r3 << std::endl;
//	  delete pb;
//  }
}


//#include "MatDataSet.h"
//#include "Allocator.h"
//using namespace Torch;
//
////@data
////0,1,0,0,1,1,1,0
////1,1,0,0,1,1,1,1
////2,1,0,0,1,1,1,0
////3,1,0,0,1,1,1,1
////4,1,1,0,1,1,1,0
////5,1,1,0,1,1,1,1
////6,1,1,0,1,1,1,0
////7,1,1,0,1,1,1,1
////@euclidean distances
////10.58,2,2,2,2,2,2.83,0
//int IOUnitTest(const string &fname, int fformat) {
//	Allocator *torchAllocator= NULL;
//	string fileName(fname);
//
//	//=================== load dataset =========================================
//	MatDataSet	*mat_data = new(torchAllocator) MatDataSet(fileName.c_str(), -1, 1, false, -1, false, fformat);
//
//	//===================  The class format ======================================  
//	//OneHotClassFormat *class_format= NULL;
//	//class_format = new(torchAllocator) OneHotClassFormat(targetDims);
//	
//	real *src_ = NULL;
//	real z;
//	real sum0=0;
//	real sum1=0;
//	assert(mat_data->n_inputs==7);	// 7 attr. 
//	assert(mat_data->n_targets==1);	// 1 target
//	assert(mat_data->n_examples==8);// 8 lines
//	assert(1);
//	for (int t = 0; t < mat_data->n_examples; t++) {
//		mat_data->setExample(t);
//		// inputs
//		src_ = mat_data->inputs->frames[0];
//		for (int j = 0; j < mat_data->n_inputs; j++) {
//			z = src_[j];
//			if(j==0) sum0 += z;
//			if(j==1) sum1 += z;
//		}
//		// targets
//		src_ = mat_data->targets->frames[0];
//		z = src_[0];
//		assert(z==(t%2));
//	}
//	assert(sum0==28);
//	assert(sum1==8);
//	return 1;
//}