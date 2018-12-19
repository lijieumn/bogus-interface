#include "RevCoulombLaw.hh"

#include <bogus/Core/Utils/NumTraits.hpp>
#include <bogus/Core/Utils/LinearSolverBase.hpp>
#include <bogus/Core/Utils/NonSmoothNewton.impl.hpp>
#include <bogus/Core/Utils/Polynomial.impl.hpp>

#include <bogus/Extra/SOC/FischerBurmeister.impl.hpp>

using namespace bogus ;

namespace argus {

// USEFUL CLASSES

//! Fischer-Burmeister function and jacobian computation, with optional change of variable
template< DenseIndexType Dimension, typename Scalar >
class RevFischerBurmeister
{

public:
  typedef bogus::LocalProblemTraits< Dimension, Scalar > Traits ;
  typedef bogus::FBBaseFunction< Dimension, Scalar > BaseFunction ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  //! Constructs an object modeling the function \f$ f : x \mapsto FB \left( mu, scaling \times x, A x + b \right) \f$
  RevFischerBurmeister(
    const Scalar mu,
    const Matrix& A,
    const Vector& b,
    const Scalar scaling )
      : m_mu( mu ), m_scaling( scaling ), m_A( A ), m_b( b )
  {}

  //! Sets \f$ fb := f(x) \f$
  void compute( const Vector& x, Vector& fb ) const ;
  //! Sets \f$ fb := f(x) \f$ and \f$ dFb\_dx := \frac {dF} {dx} \f$
  void computeJacobian( const Vector& x, Vector& fb, Matrix& dFb_dx ) const ;

private:
  Scalar m_mu ;
  Scalar m_scaling ;
  const Matrix& m_A ;
  const Vector& m_b ;

} ;

template< DenseIndexType Dimension, typename Scalar,
           local_soc_solver::Strategy Strat = local_soc_solver::RevHybrid  >
struct RevCoulombSolver
{

  typedef LocalProblemTraits< Dimension, Scalar > Traits ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  static Scalar solve(
          const typename Traits::Matrix &A,
          const typename Traits::Vector &b,
          typename Traits::Vector &x,
          const Scalar mu, const Scalar tol,
          const Scalar scaling = 1
          ) ;

} ;

// No analytic solution in the general case
template < DenseIndexType Dimension, typename Scalar >
struct AnalyticLocalRCSolver
{
	typedef LocalProblemTraits< Dimension, Scalar > Traits ;
	typedef typename Traits::Vector Vector ;
	typedef typename Traits::Matrix Matrix ;

	static bool solveOrthogonality(
	        const typename Traits::Matrix &,
	        const typename Traits::Vector &,
	        typename Traits::Vector &,
	        const Scalar
	        )
	{
		return false ;
	}
} ;


//SOLVER  LOGIC



template < unsigned Dimension, typename Scalar, local_soc_solver::Strategy Strat >
bool RevCoulombLaw< Dimension, Scalar, Strat>::solveLocal(const unsigned problemIndex,
            const typename Traits::Matrix &A,
            const typename Traits::Vector &b,
            typename Traits::Vector &xm , const Scalar scaling ) const
{
	typedef RevCoulombSolver< Traits::dimension, typename Traits::Scalar, local_soc_solver::PureNewton > LocalSolver ;
	Scalar err = LocalSolver::solve(  A, b, xm, m_mu[ problemIndex ], m_localTol, scaling ) ;
	return m_localTol > err;
}


template< DenseIndexType Dimension, typename Scalar, local_soc_solver::Strategy Strat >
Scalar RevCoulombSolver< Dimension, Scalar, Strat >::solve(
        const typename Traits::Matrix &A0,
        const typename Traits::Vector &b0,
        typename Traits::Vector &x,
        const Scalar mu, const Scalar tol, const Scalar scaling
        )
{
	const typename Traits::Matrix A = A0 * scaling ;
	const typename Traits::Vector b = b0 * scaling ;

	if(mu < 0 )
	{
		x = A.ldlt().solve( -b ) ;
		return (A*x + b).squaredNorm() ;
	}

	// see [Daviet et al 2011], Appendix B.2

	// Newton solver
	typedef RevFischerBurmeister< Dimension, Scalar > FBFunc ;
	typedef typename Traits::LUType LUType ;

	FBFunc fb( mu, A, b, 1 ) ;
	NonSmoothNewton< FBFunc > nsNewton( fb, tol )  ;

	// Sticking case -- u = 0, r \in Kmu
	if( mu * Traits::np( b ) >= Traits::tp( b ).norm() )
	{
		x.setZero() ;
		return 0 ;
	}

	double res = 0. ;

	if( Strat == local_soc_solver::Hybrid )
	{
		res = nsNewton.solve( x ) ;
		if( res < tol ) return res ;
	}

	// Continuing enumerative fallback

	// Take off -- r = 0, u_N > 0
	Vector x0 = x ;
	LUType( A ).solve( -b, x ) ;
	if ( Traits::np(x) >= 0.) {
		return 0. ;
	}
	// Sliding case

	if( NumTraits< Scalar >::isZero( mu ) )
	{
		// mu = 0 case
		Traits::np(x) = 0 ;
		Traits::tp(x) = A.template block<Dimension-1, Dimension -1>(1,1).fullPivLu().solve( -Traits::tp(b) ) ;
		return  ( A.template block<Dimension-1, Dimension -1>(1,1)*Traits::tp(x) + Traits::tp(b) ).squaredNorm() ;
	}


	// Sliding case
	if( ! AnalyticLocalRCSolver< Dimension, Scalar >::solveOrthogonality( A, b, x, mu ) )
	{
		x = x0 ;
	}

	// Refinement of final solution
	if( Strat == local_soc_solver::PureNewton || Strat == local_soc_solver::RevHybrid ) {
		res = nsNewton.solve( x ) ;
	} else if( Strat == local_soc_solver::Hybrid  ) {
		const double refinedRes = nsNewton.solve( x ) ;
		if( refinedRes <= res )
			return refinedRes ;

		//This can happen if the quartic solver returned a very bad value, like an
		// unreastically huge alpha
		x = x0 ;
	}

	return res ;
}

// FB NEWTON

template< DenseIndexType Dimension, typename Scalar >
void RevFischerBurmeister< Dimension, Scalar >::compute(
    const Vector& x, Vector& fb ) const
{
  const Vector u = m_scaling * x ;
  const Vector r = m_A * x + m_b ;
  FischerBurmeister< Dimension, Scalar, true >::compute( m_mu, r, u, fb ) ;
}

template< DenseIndexType Dimension, typename Scalar >
void RevFischerBurmeister< Dimension, Scalar >::computeJacobian(
      const Vector& x, Vector& fb, Matrix& dFb_dx ) const
{
  Vector ut = m_scaling * x ;
  Scalar s = Traits::tp( ut ).norm() ;
  Traits::np( ut ) += m_mu * s ;
  Vector r ( m_A * x + m_b ) ;


  Matrix dFb_dr ;
  BaseFunction::computeJacobian( m_mu, r, ut, fb, dFb_dr, dFb_dx ) ;

  if ( !NumTraits< Scalar >::isZero( s ) )
  {
	Traits::tc( dFb_dx).noalias() +=
	  Traits::nc( dFb_dx ) *  ( m_mu / s ) * Traits::tp( ut ).transpose() ;

  }

  dFb_dx *= m_scaling ;
  dFb_dx.noalias() += dFb_dr * m_A ;
}

// Specialization for Coulomb 3D Friction
template< typename Scalar >
struct AnalyticLocalRCSolver< 3u, Scalar>
{
	enum { Dimension = 3 } ;
	typedef LocalProblemTraits< Dimension, Scalar > Traits ;
	typedef typename Traits::Vector Vector ;
	typedef typename Traits::Matrix Matrix ;

	static bool solveOrthogonality(
	        const typename Traits::Matrix &A,
	        const typename Traits::Vector &b,
	        typename Traits::Vector &u,
	        const Scalar mu
	        )
	{
		// By Liam Toran

		// mu > 0 case
		// Calculation aids
		Scalar Lambda = A(2,2)*b(1) - A(1,2)*b(2);
		Scalar Gamma = A(1,1)*b(2) - A(1,2)*b(1);
		Scalar detAt = A(1,1)*A(2,2) - A(1,2)*A(1,2);
		Scalar trAt = A(1,1)+ A(2,2);
		Scalar Psi = A(0,1)*Lambda + A(0,2)*Gamma - b(0)*detAt;
		Scalar Omega = A(0,1)*b(1) + A(0,2)*b(2) - b(0)*trAt;
		Scalar D = - b(0);
		Scalar coeff_dom = -(mu*mu)*(Psi*Psi);
		Scalar c3_aux = -2*(mu*mu)*Psi*Omega;
		Scalar c2_aux = Lambda*Lambda + Gamma*Gamma -(mu*mu)*(2*Psi*D + Omega*Omega);
		Scalar c1_aux = 2*(Lambda*b(1) + Gamma*b(2) - ((mu*mu)*D*Omega));
		Scalar c0_aux = b(1)*b(1) + b(2)*b(2) - (mu*D)*(mu*D);
		// Coefficient of the quartic equation
		// c0 + c1*alpha + ... + c3*alpha**3 + alpha **4 = 0
		Scalar c0 = c0_aux/coeff_dom;
		Scalar c1 = c1_aux/coeff_dom;
		Scalar c2 = c2_aux/coeff_dom;
		Scalar c3 = c3_aux/coeff_dom;
		Scalar coeffs[4] = { c0 , c1 , c2, c3 } ;
		// Solving the quartic equation
		Scalar roots[4] ;
		unsigned nRoots = bogus::polynomial::getRealRoots( coeffs, roots, bogus::polynomial::StrictlyPositiveRoots ) ;


		//Getting a solution to the Coulomb Frictionnal problem.
		//We have to make sure that our root is strictly positive and that Au + b = r and mu*rn = rt.norm()
		Scalar alpha = 0. ;
		for( unsigned k = 0 ; k < nRoots ; ++k ) {
			Scalar alpha_test = roots[k];
			Scalar detG_test=  (alpha_test * alpha_test) * detAt + alpha_test * trAt + 1. ;
			Scalar rt1_test = ( alpha_test * Lambda + b(1) ) / detG_test ;
			Scalar rt2_test = ( alpha_test * Gamma + b(2) ) / detG_test ;
			Scalar u1_test = -alpha_test * rt1_test ;
			Scalar u2_test = -alpha_test * rt2_test ;
			Scalar rn_test = b(0) + A(0,1)*u1_test + A(0,2)*u2_test;
			if (rn_test>=0 && alpha > alpha_test) {
				alpha = alpha_test; // This is a solution to the Coulomb Problem.
			}
		}
		if (alpha > 0.) {
			Scalar detG =  (alpha * alpha) * detAt + alpha * trAt + 1. ;
			u[0] = 0 ;
			Scalar rt1 = ( alpha * Lambda + b(1) ) / detG ;
			Scalar rt2 = ( alpha * Gamma + b(2) ) / detG ;
			u[1] = - alpha * rt1 * b.norm() / A.norm();
			u[2] = - alpha * rt2 * b.norm() / A.norm();

			return true ;
		}

		// There is no solution to the Coulomb Problem
		return false ;
	}

} ;

} //argus

