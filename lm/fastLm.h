// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// 
// Copyright (C) 2018 Fabio Rosa
//
// This file is an adaptation from the Rcpp/Eigen example of a simple lm() file "fastLm.h"
// (Copyright (C) Bates, Eddelbuettel Francois)
// This version does not include the GESDD solver (due to F77 deps)
// This version does not include/link dependencies from s) from Rcpp, RcppEigen or R_ext
//
// RcppEigen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppEigen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// EigenUtils is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// EigenUtils is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You can check a copy of the GNU General Public License at
// <http://www.gnu.org/licenses/>.

#ifndef RCPPEIGEN_FASTLM_H
#define RCPPEIGEN_FASTLM_H

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include <tuple>

namespace lmsol
{
	using Eigen::ArrayXd;
	using Eigen::ColPivHouseholderQR;
	using Eigen::ComputeThinU;
	using Eigen::ComputeThinV;
	using Eigen::HouseholderQR;
	using Eigen::JacobiSVD;
	using Eigen::LDLT;
	using Eigen::LLT;
	using Eigen::Lower;
	using Eigen::Map;
	using Eigen::MatrixXd;
	using Eigen::SelfAdjointEigenSolver;
	using Eigen::SelfAdjointView;
	using Eigen::TriangularView;
	using Eigen::VectorXd;
	using Eigen::Upper;

	typedef MatrixXd::Index Index;
	typedef MatrixXd::Scalar Scalar;
	typedef MatrixXd::RealScalar RealScalar;
	typedef ColPivHouseholderQR<MatrixXd>::PermutationType Permutation;
	typedef Permutation::IndicesType Indices;


	class lm
	{
	protected:
		MatrixXd m_X;                     /**< model matrix */
		VectorXd m_y;                     /**< response vector */
		Index m_n;                        /**< number of rows of X */
		Index m_p;                        /**< number of columns of X */
		VectorXd m_coef;                  /**< coefficient vector */
		int m_r;                          /**< computed rank or NA_INTEGER */
		VectorXd m_fitted;                /**< vector of fitted values */
		VectorXd m_se;                    /**< standard errors  */
		RealScalar m_prescribedThreshold; /**< user specified tolerance */
		bool m_usePrescribedThreshold;

	public:
		lm(const MatrixXd&, const VectorXd&);

		ArrayXd Dplus(const ArrayXd& D);
		MatrixXd I_p() const
		{
			return MatrixXd::Identity(m_p, m_p);
		}
		MatrixXd XtX() const;

		// setThreshold and threshold are based on ColPivHouseholderQR methods
		RealScalar threshold() const;
		const VectorXd& se() const
		{
			return m_se;
		}
		const VectorXd& coef() const
		{
			return m_coef;
		}
		const VectorXd& fitted() const
		{
			return m_fitted;
		}
		int rank() const
		{
			return m_r;
		}
		lm& setThreshold(const RealScalar&);
	};

	class ColPivQR : public lm
	{
	public:
		ColPivQR(const MatrixXd&, const VectorXd&);
	};

	class Llt : public lm
	{
	public:
		Llt(const MatrixXd&, const VectorXd&);
	};

	class Ldlt : public lm
	{
	public:
		Ldlt(const MatrixXd&, const VectorXd&);
	};

	class QR : public lm
	{
	public:
		QR(const MatrixXd&, const VectorXd&);
	};

	class SVD : public lm
	{
	public:
		SVD(const MatrixXd&, const VectorXd&);
	};

	class SymmEigen : public lm
	{
	public:
		SymmEigen(const MatrixXd&, const VectorXd&);
	};

	enum { ColPivQR_t = 0, QR_t, LLT_t, LDLT_t, SVD_t, SymmEigen_t, GESDD_t };
	enum { coefficients = 0, se, rank, df_residual, residuals, s, fitted_values };
	typedef std::tuple<VectorXd, VectorXd, int, int, VectorXd, double, VectorXd> lmres;
	//        		 List::create(_["coefficients"]  = coef,
	//                            _["se"]            = se,
	//                            _["rank"]          = rank,
	//                            _["df.residual"]   = df,
	//                            _["residuals"]     = resid,
	//                            _["s"]             = s,
	//                            _["fitted.values"] = ans.fitted());

	// extern "C" SEXP fastLm(SEXP Xs, SEXP ys, SEXP types);
	lmres fastLm(const MatrixXd&, const VectorXd&, int type);
}

namespace lmsol
{

	//double NA_REAL;
	//int NA_INTEGER;
#define NA_REAL 0.0
#define NA_INTEGER 0

	using std::invalid_argument;
	using std::numeric_limits;

	lm::lm(const MatrixXd& X, const VectorXd& y)
		: m_X(X)
		, m_y(y)
		, m_n(X.rows())
		, m_p(X.cols())
		, m_coef(VectorXd::Constant(m_p, NA_REAL))
		, m_r(NA_INTEGER)
		, m_fitted(m_n)
		, m_se(VectorXd::Constant(m_p, NA_REAL))
		, m_usePrescribedThreshold(false)
	{
	}

	lm& lm::setThreshold(const RealScalar& threshold)
	{
		m_usePrescribedThreshold = true;
		m_prescribedThreshold = threshold;
		return *this;
	}

	inline ArrayXd lm::Dplus(const ArrayXd& d)
	{
		ArrayXd di(d.size());
		double comp(d.maxCoeff() * threshold());
		for (int j = 0; j < d.size(); ++j)
			di[j] = (d[j] < comp) ? 0. : 1. / d[j];
		m_r = (di != 0.).count();
		return di;
	}

	MatrixXd lm::XtX() const
	{
		return MatrixXd(m_p, m_p).setZero().selfadjointView<Lower>().rankUpdate(m_X.adjoint());
	}

	/** Returns the threshold that will be used by certain methods such as rank().
	 *
	 *  The default value comes from experimenting (see "LU precision
	 *  tuning" thread on the Eigen list) and turns out to be
	 *  identical to Higham's formula used already in LDLt.
	 *
	 *  @return The user-prescribed threshold or the default.
	 */
	RealScalar lm::threshold() const
	{
		return m_usePrescribedThreshold ? m_prescribedThreshold : numeric_limits<double>::epsilon() * m_p;
	}

	ColPivQR::ColPivQR(const MatrixXd& X, const VectorXd& y)
		: lm(X, y)
	{
		ColPivHouseholderQR<MatrixXd> PQR(X); // decompose the model matrix
		Permutation Pmat(PQR.colsPermutation());
		m_r = PQR.rank();
		if (m_r == m_p) { // full rank case
			m_coef = PQR.solve(y);
			m_fitted = X * m_coef;
			m_se = Pmat * PQR.matrixQR().topRows(m_p).triangularView<Upper>().solve(I_p()).rowwise().norm();
			return;
		}
		MatrixXd Rinv(PQR.matrixQR().topLeftCorner(m_r, m_r).triangularView<Upper>().solve(MatrixXd::Identity(m_r, m_r)));
		VectorXd effects(PQR.householderQ().adjoint() * y);
		m_coef.head(m_r) = Rinv * effects.head(m_r);
		m_coef = Pmat * m_coef;
		// create fitted values from effects
		// (can't use X*m_coef if X is rank-deficient)
		effects.tail(m_n - m_r).setZero();
		m_fitted = PQR.householderQ() * effects;
		m_se.head(m_r) = Rinv.rowwise().norm();
		m_se = Pmat * m_se;
	}

	QR::QR(const MatrixXd& X, const VectorXd& y)
		: lm(X, y)
	{
		HouseholderQR<MatrixXd> QR(X);
		m_coef = QR.solve(y);
		m_fitted = X * m_coef;
		m_se = QR.matrixQR().topRows(m_p).triangularView<Upper>().solve(I_p()).rowwise().norm();
	}

	Llt::Llt(const MatrixXd& X, const VectorXd& y)
		: lm(X, y)
	{
		LLT<MatrixXd> Ch(XtX().selfadjointView<Lower>());
		m_coef = Ch.solve(X.adjoint() * y);
		m_fitted = X * m_coef;
		m_se = Ch.matrixL().solve(I_p()).colwise().norm();
	}

	Ldlt::Ldlt(const MatrixXd& X, const VectorXd& y)
		: lm(X, y)
	{
		LDLT<MatrixXd> Ch(XtX().selfadjointView<Lower>());
		Dplus(Ch.vectorD()); // to set the rank
							 // FIXME: Check on the permutation in the LDLT and incorporate it in
		// the coefficients and the standard error computation.
		//	m_coef            = Ch.matrixL().adjoint().
		//	    solve(Dplus(D) * Ch.matrixL().solve(X.adjoint() * y));
		m_coef = Ch.solve(X.adjoint() * y);
		m_fitted = X * m_coef;
		m_se = Ch.solve(I_p()).diagonal().array().sqrt();
	}

	SVD::SVD(const MatrixXd& X, const VectorXd& y)
		: lm(X, y)
	{
		JacobiSVD<MatrixXd> UDV(X.jacobiSvd(ComputeThinU | ComputeThinV));
		MatrixXd VDi(UDV.matrixV() * Dplus(UDV.singularValues().array()).matrix().asDiagonal());
		m_coef = VDi * UDV.matrixU().adjoint() * y;
		m_fitted = X * m_coef;
		m_se = VDi.rowwise().norm();
	}

	SymmEigen::SymmEigen(const MatrixXd& X, const VectorXd& y)
		: lm(X, y)
	{
		SelfAdjointEigenSolver<MatrixXd> eig(XtX().selfadjointView<Lower>());
		MatrixXd VDi(eig.eigenvectors() * Dplus(eig.eigenvalues().array()).sqrt().matrix().asDiagonal());
		m_coef = VDi * VDi.adjoint() * X.adjoint() * y;
		m_fitted = X * m_coef;
		m_se = VDi.rowwise().norm();
	}

	static inline lm do_lm(const MatrixXd& X, const VectorXd& y, int type)
	{
		switch (type) {
		case ColPivQR_t:
			return ColPivQR(X, y);
		case QR_t:
			return QR(X, y);
		case LLT_t:
			return Llt(X, y);
		case LDLT_t:
			return Ldlt(X, y);
		case SVD_t:
			return SVD(X, y);
		case SymmEigen_t:
			return SymmEigen(X, y);
			//    case GESDD_t:
			//	return GESDD(X, y);
		}
		throw invalid_argument("invalid type");
		return ColPivQR(X, y); // -Wall
	}

	lmres fastLm(const MatrixXd& X, const VectorXd& y, int type)
	{
		Index n = X.rows();
		if (y.size() != n)
			throw invalid_argument("size mismatch");

		// Select and apply the least squares method
		lm ans(do_lm(X, y, type));

		// Copy coefficients and install names, if any
		VectorXd coef = ans.coef();
		VectorXd resid = y - ans.fitted();
		int rank = ans.rank();
		int df = (rank == NA_INTEGER) ? n - X.cols() : n - rank;
		double s = resid.norm() / std::sqrt(double(df));

		// Create the standard errors
		VectorXd se = s * ans.se();
		VectorXd fit = ans.fitted();

		lmres res = std::make_tuple(coef, se, rank, df, resid, s, fit);
		return res;
	}
}



#endif
