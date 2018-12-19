/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_IMPL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_IMPL_HPP

#include "NodalGaussSeidel.hh"

#include <bogus/Core/BlockSolvers/GaussSeidelBase.impl.hpp>
#include <bogus/Core/BlockSolvers/Coloring.impl.hpp>
#include <bogus/Core/BlockSolvers/Krylov.impl.hpp>

#ifndef BOGUS_DONT_PARALLELIZE
#include <omp.h>
#endif

#include <numeric>
#include <fstream>
#include <sstream>

#define LAMBDA_FROM_FORCES 1
#define OUTPUT 0

namespace argus
{

template < typename BlockMatrixType, typename CstMatrixType >
NodalGaussSeidel< BlockMatrixType, CstMatrixType >& NodalGaussSeidel< BlockMatrixType, CstMatrixType >::setMatrix(
        const bogus::BlockObjectBase< BlockMatrixType > & M )
{
	if( m_matrix != &M && ( m_matrix != BOGUS_NULL_PTR( const BlockObjectBase< BlockMatrixType >) ||
	                        m_coloring.size() != (std::size_t) M.rowsOfBlocks() )) {
		m_coloring.update( false, M.derived() );
	}

	m_matrix = &M ;

	updateLocalMatrices() ;

	return *this ;
}

template < typename BlockMatrixType, typename CstMatrixType >
void NodalGaussSeidel< BlockMatrixType, CstMatrixType >::updateLocalMatrices( )
{

	if( !m_matrix )
		return ;

	const Index n = m_matrix->rowsOfBlocks() ;
	m_localMatrices.resize( n ) ;

	m_force_scaling.resize( n ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( Index i = 0 ; i <  n ; ++i )
	{
		const typename BlockMatrixType::BlockPtr ptr = Base::explicitMatrix().diagonalBlockPtr( i ) ;

		if( ptr == BlockMatrixType::InvalidBlockPtr ) {
			bogus::resize( m_localMatrices[i], m_matrix->blockRows(i), m_matrix->blockCols(i) ) ;
			bogus::set_zero( m_localMatrices[i] ) ;
			m_force_scaling[i] = 1. ;
		} else {
			m_localMatrices[i] = bogus::MatrixTraits<typename BlockMatrixType::BlockType>
			        ::asConstMatrix( Base::explicitMatrix().block( ptr ) ) ;
			if (m_mass_inv.empty()) {
				m_force_scaling[i] = std::max(Scalar(1), m_localMatrices[i].rows()
				                           / (bogus::NumTraits<Scalar>::epsilon() + m_localMatrices[i].trace()) ) ;
			} else {
				m_force_scaling[i] = m_mass_inv[i];
			}
		}
	}

//	if(! m_force_scaling.empty())
//		std::cout << "Avg force scaling: " << std::accumulate(m_force_scaling.begin(), m_force_scaling.end(), 0.0) / m_force_scaling.size() << std::endl ;

	m_lambda_scaling =
	        std::accumulate(m_force_scaling.begin(), m_force_scaling.end(), 0.)
	        / (1+m_force_scaling.size()) ;

	Base::processLocalMatrices() ;
}


template < typename BlockMatrixType, typename CstMatrixType >
template < typename NSLaw,  typename RhsT, typename ResT, typename LambdaT, typename DeltaLambdaT >
void NodalGaussSeidel< BlockMatrixType, CstMatrixType >::innerLoop(
        bool parallelize, const NSLaw &law,
        const RhsT& b,
        DeltaLambdaT& r,
        std::vector< unsigned char > &skip,
        Scalar &ndxRef,
        ResT &x,
        LambdaT& lambda	) const
{
	typedef typename NSLaw::Traits LocalProblemTraits ;
	const Index dimension = Base::BlockProblemTraits::dimension ;

	ResT bl = b+lambda ;

	bogus::Segmenter< dimension, ResT, typename BlockMatrixType::Index >
	        xSegmenter( x, m_matrix->rowOffsets() ) ;
	const bogus::Segmenter< dimension, const RhsT, typename BlockMatrixType::Index >
	        bSegmenter( b, m_matrix->rowOffsets() ) ;
	bogus::Segmenter< dimension, LambdaT, typename BlockMatrixType::Index >
	        lSegmenter( lambda, m_matrix->rowOffsets() ) ;
	bogus::Segmenter< dimension, DeltaLambdaT, typename BlockMatrixType::Index >
	        rSegmenter( r, m_matrix->rowOffsets() ) ;

	const Scalar absSkipIters = std::min( m_skipIters, (unsigned) std::sqrt( (Scalar) skip.size() ) ) ;

#if OUTPUT
	std::ofstream out;
	static int frame = 0;
	if (frame == 0) {
		out.open("output_toplot.txt");
	} else {
		out.open("output_toplot.txt", std::ios::app);
	}
	Scalar local_err = 0;
#endif


#ifdef BOGUS_DONT_PARALLELIZE
	(void) parallelize ;
#else
#pragma omp parallel if ( parallelize )
	{
#endif
		typename LocalProblemTraits::Vector lb, lx, ldx, ldl ;

		for( unsigned c = 0 ; c+1 < m_coloring.colors.size() ; ++ c )
		{

			Scalar locNdxRef = ndxRef ;

#ifndef BOGUS_DONT_PARALLELIZE
#if OUTPUT
#pragma omp for reduction(+:local_err)
#else
#pragma omp for
#endif
#endif
			for( std::ptrdiff_t pi = m_coloring.colors[c] ; pi < m_coloring.colors[c+1] ; ++ pi )
			{

				const std::size_t i = m_coloring.permutation[ pi ] ;

				if( skip[i] ) {
					--skip[i] ;
					continue ;
				}

				lx = xSegmenter[ i ] ;
				lb = bSegmenter[ i ] + lSegmenter[ i ] - m_regularization(i) * lx ;
				Base::explicitMatrix().splitRowMultiply( i, x, lb ) ;
				ldx = -lx ;

				const bool ok = law.solveLocal( i, m_localMatrices[i], lb, lx, m_force_scaling[ i ] ) ;
				ldx += lx ;
				// local_err += law.eval(i, lx, m_localMatrices[i]*lx + lb);
#if OUTPUT
				local_err += ldx.norm();
#endif
				if( !ok ) { ldx *= .5 ;}
				xSegmenter[ i ] += ldx ;

#if LAMBDA_FROM_FORCES
//				ldl = -lSegmenter[i] ;
//				lSegmenter[ i ] = m_localMatrices[i]*xSegmenter[i] + lb  ;
				ldl = -rSegmenter[i];
				rSegmenter[i] = m_localMatrices[i]*xSegmenter[i] + lb  ;
				ldl += rSegmenter[i];
				m_constraintMat.derived().template colMultiply<false>( i, m_constraintStepSize*ldl, lambda ) ;
#endif

				const Scalar fs2 = m_force_scaling[i]*m_force_scaling[i] ;
				const Scalar nr2 = fs2 * rSegmenter[i].squaredNorm() ;
				const Scalar ndx2 = ldx.squaredNorm() + fs2*ldl.squaredNorm() ;

				if( ndx2 > locNdxRef ) locNdxRef = ndx2 ;

				if( ndx2 < m_skipTol * std::max( nr2, locNdxRef ) )
				{
					skip[i] = absSkipIters ;
				}
			}

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp critical
#endif
			if( locNdxRef > ndxRef ) ndxRef = locNdxRef ;
		}

#ifndef BOGUS_DONT_PARALLELIZE
	}
#endif

#if OUTPUT
	if ((frame++)%25 == 0) {
		out << local_err << std::endl;
	}
#endif
}

template < typename BlockMatrixType, typename CstMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename NodalGaussSeidel< BlockMatrixType, CstMatrixType >::Scalar
NodalGaussSeidel< BlockMatrixType, CstMatrixType >::solve( const NSLaw &law,
                                   const RhsT &b, ResT &x, ResT &r) const
{

	typename GlobalProblemTraits::DynVector lambda, s ;
	lambda.setZero(b.rows()) ;
	s.setZero(x.rows()) ;

	return solve(law, b, x, r, s, lambda) ;
}

template < typename BlockMatrixType, typename CstMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename NodalGaussSeidel< BlockMatrixType, CstMatrixType >::Scalar
NodalGaussSeidel< BlockMatrixType, CstMatrixType >::solve( const NSLaw &law,
                                   const RhsT &b, ResT &x) const
{
	ResT r(x.size());
	return solve(law, b, x, r);
}

template < typename BlockMatrixType, typename CstMatrixType >
template < typename NSLaw, typename RhsT, typename ResT,
           typename CstRhs, typename CstForce >
typename NodalGaussSeidel< BlockMatrixType, CstMatrixType >::Scalar
NodalGaussSeidel< BlockMatrixType, CstMatrixType >::solve(
        const NSLaw &law,
        const RhsT &b, ResT &x,
        const CstRhs& s, CstForce& lambda ) const
{
	ResT r(x.size());
	return solve(law, b, x, r, s, lambda);
}

template < typename BlockMatrixType, typename CstMatrixType >
template < typename NSLaw, typename RhsT, typename ResT,
           typename CstRhs, typename CstForce >
typename NodalGaussSeidel< BlockMatrixType, CstMatrixType >::Scalar
NodalGaussSeidel< BlockMatrixType, CstMatrixType >::solve(
        const NSLaw &law,
        const RhsT &b, ResT &x, ResT &r,
        const CstRhs& s, CstForce& lambda ) const
{
	(void) s;
	assert( m_matrix ) ;

	typename GlobalProblemTraits::DynVector x_best, z, delta_lambda (lambda.rows()) ;


	r = lambda + b;
	m_matrix->template multiply< false >( x, r, 1, 1 ) ;
	delta_lambda = lambda ;

#if LAMBDA_FROM_FORCES
	m_constraintMat.template multiply<false>( r, delta_lambda, 1, -1) ;
#else
	delta_lambda = m_constraintMat * ( b + (*m_matrix) * (x-s) ) ;
#endif
	lambda += delta_lambda;

	Scalar err_best = std::numeric_limits< Scalar >::max() ;
	evalAndKeepBest( law, x, r, delta_lambda, x_best, err_best, false) ;

	this->m_callback.trigger( 0, err_best ) ;

	const Index n = m_matrix->rowsOfBlocks() ;

	bogus::WithMaxThreads wmt( m_maxThreads ) ;
	const int newMaxThreads = wmt.nThreads() ;
	const bool parallelize = (newMaxThreads != 1 && n > newMaxThreads*newMaxThreads ) ;

	std::vector< unsigned char > skip( n, 0 ) ;
	Scalar ndxRef = 0 ; //Reference step size

	unsigned GSIter ;
	for( GSIter = 1 ; GSIter <= m_maxIters ; ++GSIter )
	{
		// (I-C) y =  Mx + b   (x,y) in Cmu
		// lambda = Cy, y = Mx + b + lambda, (x,y) in Cmu
		// lambda = C(Mx + b +lambda), (x, Mx + b + lambda) \in Cmu
		// (I-C) lambda = C(Mx + b),  (x, Mx + b + lambda) \in Cmu

		// lambda = C * ( b + lambda + M * (x-s) )
		// (I - C) lambda = z = C(b + M (x-s)) ,
		// (I - C_jj) lambda_j = z + sum_{i!=j} C_ji lambda_i


#if LAMBDA_FROM_FORCES
		z = lambda ;
#else
		lambda += m_constraintStepSize * delta_lambda ;
#endif

		for (int aa = 0; aa < 1; aa++) {
			innerLoop( parallelize, law, b, r, skip, ndxRef, x, lambda ) ;
		}

#if LAMBDA_FROM_FORCES
		m_constraintMat.template multiply<false>( r-z, delta_lambda, 1, 0) ;
		lambda = z + m_constraintStepSize * delta_lambda ;
#else
		delta_lambda = m_constraintMat * ( b + (*m_matrix) * (x-s) ) ;
#endif

		if( 0 == ( GSIter % m_evalEvery ) )
		{
			r = b + lambda ;
			m_matrix->template multiply< false >( x, r, 1, 1 ) ;

			Scalar err = evalAndKeepBest( law, x, r, delta_lambda, x_best, err_best, false ) ;

			this->m_callback.trigger( GSIter, err ) ;

			if( err < m_tol )
			{
				break ;
			}

			ndxRef /= m_evalEvery ;
		}

	}

	if( GSIter > m_maxIters ) x = x_best ;

	r = b + lambda + (*m_matrix) * x;

	return err_best ;

}


template < typename BlockMatrixType, typename CstMatrixType >
typename NodalGaussSeidel< BlockMatrixType, CstMatrixType >::Scalar
NodalGaussSeidel< BlockMatrixType, CstMatrixType >::evalConstraintsError(
        const typename GlobalProblemTraits::DynVector& delta_lambda
) const
{
	Scalar err = -1 ;
	if( Base::m_useInfinityNorm ) {
		err = m_lambda_scaling * delta_lambda.template lpNorm<Eigen::Infinity>() ;
	} else {
		err = m_lambda_scaling * m_lambda_scaling * delta_lambda.squaredNorm() / (delta_lambda.rows()+1) ;
	}
	// std::cout << "Constraint err: " << err << std::endl ;

	return err ;
}

template < typename BlockMatrixType, typename CstMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename NodalGaussSeidel< BlockMatrixType, CstMatrixType >::Scalar
NodalGaussSeidel< BlockMatrixType, CstMatrixType >::eval(
        const NSLaw &law,
        const ResT &y, const RhsT &x, bool originalValue ) const
{
	const Index dimension = Base::BlockProblemTraits::dimension ;

	const bogus::Segmenter< dimension, const RhsT, typename BlockMatrixType::Index >
	        xSegmenter( x, m_matrix->rowOffsets() ) ;
	const bogus::Segmenter< dimension, const ResT, typename BlockMatrixType::Index >
	        ySegmenter( y, m_matrix->rowOffsets() ) ;

	typedef typename bogus::BlockMatrixTraits< BlockMatrixType >::Index Index ;

	const Index n = m_matrix->rowsOfBlocks() ;

	Scalar err = 0. ;
	typename NSLaw::Traits::Vector lx, ly ;

	if( Base::usesInfinityNorm() )
	{

    #ifndef BOGUS_DONT_PARALLELIZE
    #pragma omp parallel private( lx, ly )
    #endif
		{
			Scalar lres = 0. ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp for
#endif
			for( Index i = 0 ; i < n ; ++ i )
			{
				lx = xSegmenter[ i ] ;
				ly = ySegmenter[ i ] * m_force_scaling[i] ;
				if (originalValue) {
					lx *= m_eval_scaling*m_A_scaling/m_b_scaling;
					ly *= m_eval_scaling/m_b_scaling;
				}
				lres = std::max( law.eval( i, lx, ly ), lres ) ;
			}

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp critical
#endif
			err = std::max( err, lres ) ;
		}

		return err ;

	} else {


#ifdef BOGUS_DONT_PARALLELIZE
		for( Index i = 0 ; i < n ; ++ i )
		{
			lx = xSegmenter[ i ] ;
			ly = ySegmenter[ i ] * m_force_scaling[i];
			if (originalValue) {
				lx *= m_eval_scaling*m_A_scaling/m_b_scaling;
				ly *= m_eval_scaling/m_b_scaling;
			}
			lres = law.eval( i, lx, ly ) ;
			err += lres ;
		}
#else
		std::vector< Scalar > lerr ( omp_get_max_threads(), 0 ) ;

        #pragma omp parallel private(lx, ly)
		{
			const int tid = omp_get_thread_num() ;
#pragma omp for
			for( Index i = 0 ; i < n ; ++ i )
			{
				lx = xSegmenter[ i ] ;
				ly = ySegmenter[ i ] * m_force_scaling[i];
				lerr[tid] += law.eval( i, lx, ly ) ;
			}

#pragma omp single
			{
				const int num_threads = omp_get_num_threads() ;
				for( int i = 0 ; i < num_threads ; ++i ) {
					err += lerr[i] ;
				}
			}
		}
#endif
		return err / ( 1 + n );

	}

}

template < typename BlockMatrixType, typename CstMatrixType >
template < typename NSLaw, typename ResT >
typename NodalGaussSeidel< BlockMatrixType, CstMatrixType >::Scalar
NodalGaussSeidel< BlockMatrixType, CstMatrixType >::evalAndKeepBest(
        const NSLaw &law, const ResT &x,
        const typename GlobalProblemTraits::DynVector& y,
        const typename GlobalProblemTraits::DynVector& delta_lambda,
        typename GlobalProblemTraits::DynVector& x_best, Scalar &err_best, bool originalValue ) const
{
#if OUTPUT
	std::ofstream out;
	if (originalValue) {
		std::ostringstream os;
		os << "output_result" << m_A_scaling << "_" << m_b_scaling << ".txt";
		static int iter = 0;
		if (iter++ == 0) {
			out.open(os.str());
		} else {
			out.open(os.str(), std::ios::app);
		}
	}
#endif
	Scalar err = eval( law, y, x, originalValue ) ;
#if OUTPUT
	if (originalValue) {
		out << err << "\t";
	}
#endif

//	std::cout << "Contact err: " << err << std::endl ;
	Scalar const_err;
	if (originalValue) {
		const_err = evalConstraintsError(delta_lambda*m_eval_scaling/m_b_scaling ) ;
	} else {
		const_err = evalConstraintsError(delta_lambda ) ;
	}
#if OUTPUT
	if (originalValue) {
		out << const_err << std::endl;
	}
#endif
	err += const_err;

	if( err < err_best )
	{
		x_best = x ;
		err_best = err ;
	}

	return err ;
}

//! Same as solveCadoux, with r = W*u +b
/*! \warning requires mu > 0 */
template< unsigned Dimension, typename WType, typename BlockMatrixType, typename CstMatrixType >
static double solveCadouxVelWithConstraints(
        const WType& W,
        const typename WType::Scalar* b, const typename WType::Scalar* mu,
        NodalGaussSeidel< BlockMatrixType, CstMatrixType > &minimizer,
        typename WType::Scalar* u, const unsigned cadouxIterations,
        const bogus::Signal<unsigned, typename WType::Scalar> *callback = BOGUS_NULL_PTR(const void),
        const typename WType::Scalar tolTighten = 1.e-1
        )
{
	// Wu + b = r
	// u* = u + s n
	// Wu* + b - W(s n) = r

	typedef typename WType::Scalar Scalar ;
	const std::ptrdiff_t n = W.rowsOfBlocks() ;

	bogus::SOCLaw< Dimension, Scalar, false > socLaw	  ( n, mu ) ;

	Eigen::Map< Eigen::VectorXd > u_map ( u, W.rows() ) ;
	Eigen::Map< const Eigen::VectorXd > b_map ( b, W.rows() ) ;

	Eigen::VectorXd s( W.rows() ), Wsb( W.rows() ), ustar( u_map ), r( W.rows() ), lambda ;
	lambda.setZero( b_map.rows() ) ;

	double res = -1 ;
	const double tol = minimizer.tol() ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( std::ptrdiff_t i = 0 ; i < n ; ++i )
	{
		s[ Dimension*i ] = u_map.segment< Dimension-1 >( Dimension*i+1 ).norm() * std::max(0., 1./mu[i]) ;
		s.segment< Dimension-1  >( Dimension*i+1 ).setZero() ;
	}
	ustar = u_map + s ;

	//Evaluate intial error
	r = W * u_map + b_map + lambda ;
	res = minimizer.eval( socLaw, r, ustar ) ;
	if( callback ) callback->trigger( 0, res ) ;

	for( unsigned cdxIter = 0 ; cdxIter < cadouxIterations ; ++cdxIter )
	{
		minimizer.setTol( tolTighten * std::max( tol, std::min( res, 1. ) ) ) ;

		Wsb = b_map - W * s  ;

		minimizer.solve( socLaw, Wsb, ustar, s, lambda ) ;

		u_map = ustar - s ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( std::ptrdiff_t i = 0 ; i < n ; ++i )
		{
			s[ Dimension*i ] = ustar.segment< Dimension-1 >( Dimension*i+1 ).norm() * std::max(0., 1./mu[i]) ;
			s.segment< Dimension-1  >( Dimension*i+1 ).setZero() ;
		}

		ustar = u_map + s ;

		// Evaluate current error

		r = W * u_map + b_map + lambda ;
		res = minimizer.eval( socLaw, r, ustar ) ;

		if( callback ) callback->trigger( cdxIter+1, res ) ;
		if( cdxIter > 0 && res < tol ) break ;

	}

	minimizer.setTol( tol ) ;

	return res ;
}


} //namespace bogus


#endif
