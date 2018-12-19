/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef ARGUS_NODAL_GAUSS_SEIDEL_HPP
#define ARGUS_NODAL_GAUSS_SEIDEL_HPP

#include <bogus/Core/BlockSolvers/GaussSeidelBase.hpp>
#include <bogus/Core/BlockSolvers/Coloring.hpp>

namespace argus
{

typedef bogus::Coloring Coloring ;

//! Projected Gauss-Seidel iterative solver.
/*!
   Works by taking into account only one block-row of the system at a time, and iterating
   several times over the whole set of rows several times until convergence has been achieved.

   Each inner iteration of the algorithm will try to solve the local problem
	  \f[
		\left\{
		  \begin{array}{rcl}
			y_i^{k+1} &=& M_{i,i} x_i^{k+1}  + b_i^{k} \\
			&s.t.& law (x^{k+1},y^{k+1})
		  \end{array}
		\right.
	  \f]
	where \b k is the current global iteration, \b i the current row
	and \f[ b_i^{k} := b_i + \sum_{ j < i }{ M_{i,j}x_j^{k+1} } +  \sum_{ j > i }{ M_{i,j}x_j^{k} } \f]

   See also solve() and \cite JAJ98.
  */
template < typename BlockMatrixType, typename CstMatrixType >
class NodalGaussSeidel : public bogus::GaussSeidelBase< NodalGaussSeidel<BlockMatrixType, CstMatrixType>, BlockMatrixType >
{
public:
	typedef bogus::GaussSeidelBase< NodalGaussSeidel, BlockMatrixType > Base ;

	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	//! Constructor with the system matrix
	NodalGaussSeidel(
	        const bogus::BlockObjectBase< BlockMatrixType > & matrix,
	        const bogus::BlockObjectBase< CstMatrixType > & cstMat,
	        const std::vector<double> m_inv
	        ) : Base(), m_lambda_scaling(1),
	            m_constraintStepSize( .2 ), m_constraintIterations( 2 ),
	            m_constraintMat( cstMat ), m_A_scaling(1), m_b_scaling(1),
	            m_eval_scaling(1), m_mass_inv(m_inv)
	{  setMatrix( matrix ) ; }

	//! Sets the system matrix and initializes internal structures
	NodalGaussSeidel& setMatrix( const bogus::BlockObjectBase< BlockMatrixType > & matrix ) ;

	//! Finds an approximate solution for a constrained linear problem
	/*!
	  Stops when the residual computed in eval() is below \ref m_tol, of the number
	  of iterations exceed \ref m_maxIters

	  Implements Algorithm 1. from \cite DBB11 to solve
	   \f[
		\left\{
		  \begin{array}{rcl}
			y &=& M x + b \\
			&s.t.& law (x,y)
		  \end{array}
		\right.
	  \f]
	  \param law The (non-smooth) law that should define:
		- An error function for the local problem
		- A local solver for each row of the system ( e.g. 1 contact solver )
		\sa SOCLaw
	  \param b the const part of the right hand side
	  \param x the unknown. Can be warm-started
	  */
	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x) const ;

	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x, ResT &r) const ;

	template < typename NSLaw, typename RhsT, typename ResT,
	           typename CstRhs, typename CstForce  >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x,
	              const CstRhs& s, CstForce& lambda) const ;

	template < typename NSLaw, typename RhsT, typename ResT,
	           typename CstRhs, typename CstForce  >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x, ResT &r,
	              const CstRhs& s, CstForce& lambda) const ;
	//! Access to the current Coloring. Will be reset whenever the matrix is changed.
	/*! Determiniticy is achieved through the mean of contact coloring ;
	  contacts that do not interact directly together can chare the same color,
	  and all contacts within a given color can be solver in parallel */
	Coloring& coloring( ) { return m_coloring ; }


	void setConstraintStepSize( float css )
	{
		m_constraintStepSize = css ;
	}

	Scalar constraintStepSize() const
	{
		return m_constraintStepSize ;
	}

	void setConstraintIterations( int ci )
	{
		m_constraintIterations = ci ;
	}

	int constraintIterations() const
	{
		return m_constraintIterations ;
	}

	void setAScaling( double A_scaling) {
		m_A_scaling = A_scaling;
	}

	void setBScaling( double b_scaling) {
		m_b_scaling = b_scaling;
	}

	void setEvalScaling( double eval_scaling) {
		m_eval_scaling = eval_scaling;
	}

	void setMassInv( const std::vector<double> &m_inv) {
		m_mass_inv = m_inv;
	}

	const bogus::BlockObjectBase< CstMatrixType > & constraintMat() const {
		return m_constraintMat ;
	}

	using Base::solve ;


	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar eval(
	        const NSLaw &law,
	        const ResT &y, const RhsT &x, bool originalValue = false ) const;

	template < typename NSLaw, typename ResT >
	Scalar evalAndKeepBest(
	    const NSLaw &law, const ResT &x,
	    const typename GlobalProblemTraits::DynVector& y,
	    const typename GlobalProblemTraits::DynVector& delta_lambda,
	    typename GlobalProblemTraits::DynVector& x_best, Scalar &err_best,
	    bool originalValue = false ) const ;

	Scalar evalConstraintsError(
	        const typename GlobalProblemTraits::DynVector& delta_lambda
	        ) const ;


protected:
	void updateLocalMatrices() ;

	template < typename NSLaw,  typename RhsT, typename ResT, typename LambdaT, typename DeltaLambdaT >
	void innerLoop (
	    bool parallelize, const NSLaw &law, const RhsT& b, DeltaLambdaT& r,
	    std::vector< unsigned char > &skip, Scalar &ndxRef,
	    ResT &x, LambdaT& lambda ) const ;

	typedef typename Base::Index Index ;

	using Base::m_matrix ;
	using Base::m_maxIters ;
	using Base::m_tol ;
	using Base::m_scaling ;
	using Base::m_maxThreads ;
	using Base::m_evalEvery ;
	using Base::m_skipTol ;
	using Base::m_skipIters ;
	using Base::m_localMatrices ;
	using Base::m_regularization ;

	std::vector< Scalar > m_force_scaling ;
	Scalar m_lambda_scaling ;

	Scalar m_constraintStepSize ;
	int    m_constraintIterations ;
	const bogus::BlockObjectBase< CstMatrixType > & m_constraintMat ;

	Scalar m_A_scaling;
	Scalar m_b_scaling;
	Scalar m_eval_scaling;
	std::vector< Scalar > m_mass_inv;

	Coloring m_coloring ;

	typedef typename Base::BlockProblemTraits::LUType LocalFacType ;
} ;

} //namespace bogus


#endif
