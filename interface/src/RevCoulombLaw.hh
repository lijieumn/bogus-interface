#ifndef ARGUS_REV_COULOMB_LAW_HPP
#define ARGUS_REV_COULOMB_LAW_HPP

#include <bogus/Extra/SecondOrder.fwd.hpp>
#include <bogus/Core/Eigen/EigenProblemTraits.hpp>
#include <bogus/Extra/SOC/FischerBurmeister.hpp>
#include <bogus/Core/Utils/NumTraits.hpp>

//Same as bogus::SOCLaw, but with reverted roles for un and r

namespace argus
{

typedef bogus::DenseIndexType DenseIndexType ;


//! Non-smooth laws based on Second Order Cone complementarity. To be used within as the first argument to GaussSeidel::solve().
/*!
	\tparam Dimension the dimension of the local problem. Specializations exist form dimension 2 and 3.
	\tparam Scalar the scalar type
	\tparam Strat local_soc_solver::Strategy for solving the local problems. Unavailable for dimensions other than 2 and 3.
  */
template < unsigned Dimension, typename Scalar,
           bogus::local_soc_solver::Strategy Strat >
class RevCoulombLaw
{
public:

	typedef bogus::LocalProblemTraits< Dimension, Scalar > Traits ;
	enum{ dimension = Dimension } ;

	//! Constructor
	/*!
	  \param n the size of the global problem ( number of contacts )
	  \param mu array containing the apertures of each second order cone ( friction coefficients )
	  */
	RevCoulombLaw( const unsigned n, const double * mu )
	    : m_mu(mu), m_n(n), m_localTol( std::pow( bogus::NumTraits< Scalar >::epsilon(), .75 ) )
	{}

	void setLocalTolerance( Scalar tol )
	{
		m_localTol = tol ;
	}

	//! \return \f$ \vert fb( mu, r, u ) \vert^2_2 \f$, where fb is the SOC Fischer-Burmeister function
	Scalar eval( const unsigned problemIndex,
	             const typename Traits::Vector &u,
	             const typename Traits::Vector &r ) const
	{
		typedef bogus::FischerBurmeister< Traits::dimension, typename Traits::Scalar, true > FBFunction ;

		if(m_mu[problemIndex] < 0) return r.squaredNorm() ;

		typename Traits::Vector fb( u.rows() ) ;
		FBFunction::compute( m_mu[problemIndex], r, u, fb ) ;

		return fb.squaredNorm() ;
	}

	//! Solves the local problem
	/*!
	  \f[
		\left\{
		  \begin{array}{rcl}
			y &=& \left( A (x + s(x)) + b \right ) \\
			K_{ \frac 1 \mu } \ni y & \perp & x \in K_{ \mu }
		  \end{array}
		\right.
	  \f]
	  where \f$ \mu \f$ is \c m_mu[\p problemIndex] and \c s is the De Saxce change of variable.

	  \param scaling Used as a scaling factor for \p x when calculating the error function
	*/
	bool solveLocal(
	        const unsigned problemIndex,
	        const typename Traits::Matrix &A,
	        const typename Traits::Vector &b,
	        typename Traits::Vector &x,
	        const Scalar scaling
	        ) const ;

private:

	const double * m_mu ;
	const unsigned m_n ;
	Scalar m_localTol ;

} ;

typedef RevCoulombLaw< 3u, double, bogus::local_soc_solver::RevHybrid > RevCoulomb3D ;



} //argus

#endif
