#include "MCPLaw.hh"

#include <bogus/Core/Utils/NumTraits.hpp>

using namespace bogus ;

namespace argus {

// USEFUL CLASSES


bool MCPLaw::solveLocal(const unsigned problemIndex,
            const typename Traits::Matrix &A,
            const typename Traits::Vector &b,
            typename Traits::Vector &x, const Scalar ) const
{
	static constexpr double eps = 1.e-16;
	static constexpr int nLocalIters = 3 ; // as per paper

	for(int i = 0 ; i<nLocalIters ; ++i)
	{
		const Scalar uN = (-b[0] - A(0,1)*x[1] - A(0,2)*x[2]) ;
		x[0] = std::max( 0.0, uN / std::max(eps, A(0,0)) ) ;

		const Scalar muN = x[0]*m_mu[problemIndex];

		const Scalar uT1 = (-b[1] - A(1,0)*x[0] - A(1,2)*x[2]) ;
		x[1] = uT1 / std::max(eps, A(1,1)) ;
		if( x[1] > muN ) x[1] = muN ;
		else if( x[1] < -muN ) x[1] = -muN ;

		const Scalar uT2 = (-b[2] - A(2,0)*x[0] - A(2,1)*x[1]) ;
		x[2] = uT2 / std::max(eps, A(2,2)) ;
		if( x[2] > muN ) x[2] = muN ;
		else if( x[2] < -muN ) x[2] = -muN ;
	}

	return true ;
}

MCPLaw::Scalar MCPLaw::eval( const unsigned problemIndex,
                     const typename Traits::Vector &x,
                     const typename Traits::Vector &y ) const
{
	Traits::Vector err ;
	err[0] = x[0] + y[0] - std::sqrt( x[0]*x[0] + y[0]*y[0] ) ;

	// Using AC for tangent components...
	const Scalar muN = std::max(0.0, x[0])*m_mu[problemIndex];

	err[1] = std::max(-muN, std::min(muN, x[1] - y[1])) - x[1];
	err[2] = std::max(-muN, std::min(muN, x[2] - y[2])) - x[2];

	return err.squaredNorm() ;
}

} //argus

