#include <assert.h>

#define a 0
#define b 1
#define c 2
#define d 3

// [[Rcpp::export]]
double S_Jaccard(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = mat[a] / (mat[a] + mat[b] + mat[c]);
    return S ;
}

// [[Rcpp::export]]
double S_Dice(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = (2*mat[a]) / (2*mat[a] + mat[b] + mat[c]);
    return S ;
}

// [[Rcpp::export]]
double S_Czekanowski(const NumericVector mat) {
    return S_Dice(mat) ;
}

// [[Rcpp::export]]
double S_3wJaccard(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = (3*mat[a]) / (3*mat[a] + mat[b] + mat[c]);
    return S ;
}

// [[Rcpp::export]]
double S_NeiLi(const NumericVector mat) {
    return S_Dice(mat) ;
}

// [[Rcpp::export]]
double S_SokalSneath1(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = mat[a] / (mat[a] + 2*(mat[b] + mat[c]));
    return S ;
}

// [[Rcpp::export]]
double S_SokalMichener(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = (mat[a]+mat[d]) / (mat[a] + mat[b] + mat[c] + mat[d]);
    return S ;
}

// [[Rcpp::export]]
double S_SokalSneath2(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double St = 2*(mat[a]+mat[d]);
    double S = St / (St + mat[b] + mat[c]);
    return S ;
}

// [[Rcpp::export]]
double S_RogerTanimoto(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = (mat[a]+mat[d]) / (mat[a] + mat[d] + 2*(mat[c] + mat[b]));
    return S ;
}

// [[Rcpp::export]]
double S_Faith(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = (mat[a]+0.5*mat[d]) / (mat[a] + mat[b] + mat[c] + mat[d]);
    return S ;
}

// [[Rcpp::export]]
double S_GowerLegendre(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = (mat[a]+mat[d]) / (mat[a] + 0.5*(mat[b] + mat[c]) + mat[d]);
    return S ;
}

// [[Rcpp::export]]
double S_Intersection(const NumericVector mat) {
   
    assert(mat.size() == 4);
    return mat[a] ;
}

// [[Rcpp::export]]
double S_InnerProd(const NumericVector mat) {
   
    assert(mat.size() == 4);
    return mat[a]+mat[d] ;
}

// [[Rcpp::export]]
double S_RusselRao(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = mat[a] / (mat[a] + mat[b] + mat[c] + mat[d]);
    return S ;
}

// [[Rcpp::export]]
double S_Hamming(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = (mat[b] + mat[c]);
    return S ;
}

// [[Rcpp::export]]
double S_Euclid(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S = std::sqrt(mat[b] + mat[c]);
    return S ;
}

// [[Rcpp::export]]
double S_SquaredEuclid(const NumericVector mat) {
   
    assert(mat.size() == 4);
    double S1 = (mat[b] + mat[c]);
    return std::sqrt(S1*S1);
}