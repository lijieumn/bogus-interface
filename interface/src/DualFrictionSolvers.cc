#include "DualFrictionSolvers.hh"

#include "ClothFrictionData.hh"
#include "SolverOptions.hh"
#include "SolverStats.hh"
#include "MCPLaw.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/BlockSolvers/ProjectedGradient.impl.hpp>
#include <bogus/Core/BlockSolvers/GaussSeidel.impl.hpp>
#include <bogus/Extra/SecondOrder.impl.hpp>

#include <bogus/Interfaces/FrictionProblem.hpp>
#include <bogus/Interfaces/Cadoux.hpp>

#include <Eigen/Sparse>
#include <Eigen/LU>

namespace argus {

struct DualFrictionSolver::MFactorization
{
	typedef Eigen::SparseMatrix<double> CSR ;
	typedef bogus::SparseLU<double>     CSR_LU ;
	typedef bogus::SparseBlockMatrix< CSR_LU > SBM ;

	SBM sbm ;

	MFactorization( const CSR& csr )
	{
		std::vector< unsigned > rpb (1, (unsigned)  csr.rows()) ;
		std::vector< unsigned > cpb (1, (unsigned)  csr.cols()) ;

		sbm.setRows( rpb );
		sbm.setCols( cpb );

		sbm.insertBack(0,0).compute( csr ) ;
		sbm.finalize();

		if( sbm.block(0).factorization().info() != Eigen::Success ) {
			std::cerr << "ClothFriction: LU factorization failed" << std::endl ;
		}
	}
};

DualFrictionSolver::DualFrictionSolver(
                const ClothFrictionData& data,
                const HType &H,
                const DynVec &w
                )
        : m_data(data), m_H(H), m_w(w), m_fac( 0 )
{ factorize() ; }

DualFrictionSolver::~DualFrictionSolver()
{
	delete m_fac ;
}

void DualFrictionSolver::factorize()
{
	MFactorization::CSR csr ;
	bogus::convert( m_data.M, csr ) ;
	m_fac = new MFactorization( csr ) ;
}

double DualFrictionSolver::solve(const SolverOptions &options, DynVec &v, DynVec& r, SolverStats &stats) const
{
	typedef MFactorization::SBM MinvType ;
	typedef bogus::Product< HType, MinvType > HMinv ;
	typedef bogus::Product< HMinv, bogus::Transpose<HType> > WExpr ;

	HMinv HMi = m_H * m_fac->sbm ;
	WExpr W   = HMi * m_H.transpose() ;

	DynVec b = m_w + HMi * m_data.f ;

	r = DynVec::Zero( W.rows() ) ; // TOTO warm-start
	double res = -1 ;

	if( options.algorithm == SolverOptions::DualProjectedGradient )	{
		// Cadoux fixed-point + marix-free projected gradient solver

		bogus::ProjectedGradient< WExpr > pg ( W ) ;
		pg.setTol( options.tolerance );
		pg.useInfinityNorm( options.useInfinityNorm );
		pg.setMaxIters( options.maxIterations );
		pg.setLineSearchIterations( options.lineSearchIterations );

		if( options.projectedGradientVariant < 0 ) {
			pg.setDefaultVariant( bogus::projected_gradient::SPG );
		} else {
			pg.setDefaultVariant( (bogus::projected_gradient::Variant) options.projectedGradientVariant );
		}

		bogus::Signal<unsigned, double> callback ;
		callback.connect( stats, &SolverStats::ack );
		res = bogus::solveCadoux<3>( W, b.data(), m_data.mu.data(), pg,
		                              r.data(), options.maxOuterIterations, &callback ) ;
	} else {

		// Gauss-Seidel w/ explicit Delassus operator
		typedef bogus::DualFrictionProblem<3>::WType WType ;
		WType W ;

		// Building Delassus matrix -- highly inefficient !
		// Converting to CSR, doing the matrix mult in CSR, then back to BSR
		{
			Eigen::SparseMatrix<double> Hsnb;
			bogus::convert(m_H, Hsnb);
			Eigen::SparseMatrix<double> HsnbT = Hsnb.transpose();

			Eigen::SparseMatrix<double> MiHt = m_fac->sbm.block(0).factorization().solve(HsnbT);
			Eigen::SparseMatrix<double, Eigen::RowMajor > Wsnb = Hsnb * MiHt ;
			bogus::convert(Wsnb, W);
			W.prune(1.e-16) ;
			W.cacheTranspose();
		}

		bogus::GaussSeidel<WType> gs(W) ;
		gs.setTol( options.tolerance );
		gs.useInfinityNorm( options.useInfinityNorm );
		gs.setMaxIters( options.maxIterations );
		gs.callback().connect( stats, &SolverStats::ack );

		bogus::Coulomb3D law( m_data.mu.rows(), m_data.mu.data() ) ;
		res = gs.solve( law, b, r, false ) ;
	}

	v = m_fac->sbm * ( m_data.f + m_H.transpose() * r ) ;
	return res ;
}

ICAFrictionSolver::ICAFrictionSolver(
                const ClothFrictionData& data,
                const HType &H,
                const DynVec &w
                )
        : m_data(data), m_H(H), m_w(w)
{ }

double ICAFrictionSolver::solve(const SolverOptions &options, DynVec &v, DynVec &r, SolverStats &stats) const
{
	const int m =  m_data.M.rowsOfBlocks() ;

	DType D, Dinv ;
	D.setRows( m );
	Dinv.setRows( m );

	D.setIdentity() ;
	Dinv.setIdentity() ;

	for( int i = 0 ; i < m ; ++i ) {
		D.block(i)  = m_data.M.diagonal(i) ;
		Dinv.block(i) = m_data.M.diagonal(i).inverse() ;
	}


	typedef bogus::Product< HType, DType > HDinv ;
	typedef bogus::DualFrictionProblem<3>::WType WType ;

	HDinv HDi = m_H * Dinv ;
	WType W   = HDi * m_H.transpose() ;


	r = DynVec::Zero( W.rows() ) ; // TODO warm-start
	double res = -1 ;

	W.prune(1.e-16) ;
	W.cacheTranspose();

	bogus::GaussSeidel<WType> gs(W) ;
	gs.setTol( options.tolerance );
	gs.useInfinityNorm( options.useInfinityNorm );
	gs.setMaxIters( options.maxIterations );
	gs.callback().connect( stats, &SolverStats::ack );

	bogus::Coulomb3D law( m_data.mu.rows(), m_data.mu.data() ) ;
	argus::MCPLaw facetedLaw( m_data.mu.rows(), m_data.mu.data() ) ;

	// delta_v ~ M^-1 H' r
	DynVec delta_v  = DynVec::Zero( D.rows() ), delta_v_prev ;
	DynVec v0 = v ;

	for( unsigned k = 0 ; k < options.maxOuterIterations ; ++k )
	{
		// b = Hv^0 + w
		// u = Hdel_v +b = HD-1Hr - HD-1*M*dv + Hdv + b = Hdv + b
		DynVec b = m_w + m_H*v0 + HDi * ( (D - m_data.M) * delta_v )  ;

		if(options.faceted)
		{
			res = gs.solve( facetedLaw, b, r, false ) ;
		} else {
			res = gs.solve( law, b, r, false ) ;
		}

		const DynVec Hr = m_H.transpose() * r ;
		for( int i = 0 ; i < m ; ++i ) {
			Eigen::Vector3d c = - Hr.segment< 3 >( 3*i ) ;
			m_data.M.splitRowMultiply( i, delta_v, c ) ;
			delta_v.segment<3>( 3*i ) = -Dinv.block(i) * c ;
		}

		v = v0 + delta_v;

		res += ( m_data.M * delta_v - Hr ).squaredNorm() / (1 + m ) ;

		if( k > 0 ) {
			res += ( delta_v - delta_v_prev ).squaredNorm() / (1 + m ) ;
			std::cout << "ICA " << k << " => " << res << std::endl ;
			if(res < options.tolerance)
				break ;
		}
		delta_v_prev = delta_v ;
	}

	return res ;
}

}
