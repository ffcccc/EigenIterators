#include <assert.h>
#include <Eigen/Dense>

using namespace Eigen;

#define a 0
#define b 1
#define c 2
#define d 3

template<class _Tp>
inline _Tp S_Jaccard(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = mat[a] / (mat[a] + mat[b] + mat[c]);
    return S ;
}

template<class _Tp>
inline _Tp S_Dice(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = (2*mat[a]) / (2*mat[a] + mat[b] + mat[c]);
    return S ;
}

template<class _Tp>
inline _Tp S_Czekanowski(const Eigen::Array<_Tp, 4, 1> &mat) {
    return S_Dice(mat) ;
}

template<class _Tp>
inline _Tp S_3wJaccard(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = (3*mat[a]) / (3*mat[a] + mat[b] + mat[c]);
    return S ;
}

template<class _Tp>
inline _Tp S_NeiLi(const Eigen::Array<_Tp, 4, 1> &mat) {
    return S_Dice(mat) ;
}

template<class _Tp>
inline _Tp S_SokalSneath1(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = mat[a] / (mat[a] + 2*(mat[b] + mat[c]));
    return S ;
}

template<class _Tp>
inline _Tp S_SokalMichener(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = (mat[a]+mat[d]) / (mat[a] + mat[b] + mat[c] + mat[d]);
    return S ;
}

template<class _Tp>
inline _Tp S_SokalSneath2(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp St = 2*(mat[a]+mat[d]);
    _Tp S = St / (St + mat[b] + mat[c]);
    return S ;
}

template<class _Tp>
inline _Tp S_RogerTanimoto(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = (mat[a]+mat[d]) / (mat[a] + mat[d] + 2*(mat[c] + mat[b]));
    return S ;
}

template<class _Tp>
inline _Tp S_Faith(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = (mat[a]+0.5*mat[d]) / (mat[a] + mat[b] + mat[c] + mat[d]);
    return S ;
}

template<class _Tp>
inline _Tp S_GowerLegendre(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = (mat[a]+mat[d]) / (mat[a] + 0.5*(mat[b] + mat[c]) + mat[d]);
    return S ;
}

template<class _Tp>
inline _Tp S_Intersection(const Eigen::Array<_Tp, 4, 1> &mat) {
    return mat[a] ;
}

template<class _Tp>
inline _Tp S_InnerProd(const Eigen::Array<_Tp, 4, 1> &mat) {
    return mat[a]+mat[d] ;
}

template<class _Tp>
inline _Tp S_RusselRao(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = mat[a] / (mat[a] + mat[b] + mat[c] + mat[d]);
    return S ;
}

template<class _Tp>
inline _Tp S_Hamming(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = (mat[b] + mat[c]);
    return S ;
}

template<class _Tp>
inline _Tp S_Euclid(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S = std::sqrt(mat[b] + mat[c]);
    return S ;
}

template<class _Tp>
inline _Tp S_SquaredEuclid(const Eigen::Array<_Tp, 4, 1> &mat) {
    _Tp S1 = (mat[b] + mat[c]);
    return std::sqrt(S1*S1);
}