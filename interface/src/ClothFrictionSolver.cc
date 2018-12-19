#include "ClothFrictionSolver.hh"

#include "ClothFrictionData.hh"

#include "NodalContact.hh"
#include "NodalSelfContact.hh"
#include "TotalContact.hh"
#include "DualFrictionSolvers.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/BlockSolvers/Krylov.impl.hpp>
#include <bogus/Core/BlockSolvers/ADMM.impl.hpp>

#include <bogus/Extra/SecondOrder.impl.hpp>

namespace argus {


//! Solves the friction problem using the Alternating Minimization Algorithm
template <typename HType>
static double solveAMA(
        const ClothFrictionData& data,
        const HType & H, const DynVec &w,
        const SolverOptions &options, DynVec &v, SolverStats &stats )
{
	   bogus::DualAMA< HType > dama( H ) ;

	   dama.setDefaultVariant( bogus::admm::Standard ) ;
	   dama.setLineSearchIterations( options.lineSearchIterations );
	   dama.setTol( options.tolerance );
	   dama.setMaxIters( options.maxIterations );

	   dama.setFpStepSize(options.amaFpStepSize);
	   dama.setProjStepSize(options.amaProjStepSize);

	   bogus::Coulomb3D law( data.mu.rows(), data.mu.data() ) ;

	   // Forces - could be warmstarted ?
	   DynVec r = DynVec::Zero( H.rows() ) ;

	   //Preconditioner for striffness matrix
	   bogus::DiagonalPreconditioner< bogus::BlockObjectBase< typename ClothFrictionData::StiffnessMatrixType > > P ;
	   P.setMatrix( data.M ) ;

	   // No linear constraints
	   bogus::Zero<double> B(0,0) ;
	   Eigen::Matrix<double,0,1> k,p;

	   // Convergence info: *very* verbose !
	   dama.callback().connect( stats, &SolverStats::ack );

	   return dama.template solveWithLinearConstraints<bogus::admm::Standard>
	           ( law, data.M, B, H, P, DynVec(-data.f), k, w, v, p, r ) ;
}

double ClothFrictionSolver::solve(const SolverOptions &options, DynVec &v, SolverStats &stats ) const
{
	DynVec r;
	return solve(options, v, r, stats);
}


double ClothFrictionSolver::solve(const SolverOptions &options, DynVec &v, DynVec &r, SolverStats &stats ) const
{
	linearSolve(v);

	return contactSolve(options, v, r, stats);
}

struct AutoVelocityDuplicater {

	AutoVelocityDuplicater( const ClothFrictionData& data, DynVec & v )
	    : m_data(data), m_vel(v), m_swapped ( v.rows() != data.M.rows() )
	{
		if( m_swapped ){
			m_data.scatterVelocities( m_vel, m_orig );
			m_orig.swap(m_vel) ;
		}
	}

	~AutoVelocityDuplicater( ) {
		if( m_swapped ) {
			m_data.gatherVelocities( m_vel, m_orig );
//			std::cout << "Duplicated velocities" << std::endl ;
//			std::cout << m_vel.transpose() << std::endl ;
		}
	}

	const ClothFrictionData& m_data ;
	DynVec& m_vel ;

	const bool m_swapped ;
	DynVec m_orig ;

};

void ClothFrictionSolver::linearSolve(DynVec& v) const 
{
	AutoVelocityDuplicater( m_data, v);

	 // Compute solution of unconstrained linear system using CG
	 bogus::Krylov< ClothFrictionData::StiffnessMatrixType > kry( m_data.M ) ;
	 kry.setMaxIters( 2 * m_data.M.rows() );

	 double cgRes = kry.asCG().parallelizeRhs(false).solve( m_data.f, v ) ;
	 if(std::isnan(cgRes) || cgRes > 1.e-12) {

#ifndef SILENCE_ARGUS
	   std::cout << "CG failed, trying GMRES " << std::endl ;
#endif
	    cgRes = kry.asGMRES().setRestart(10).solve( m_data.f, v ) ;
	 }

#ifndef SILENCE_ARGUS
	 std::cout << "ClothFriction linear solve: res = " << cgRes << std::endl;
#endif
}

double ClothFrictionSolver::contactSolve(const SolverOptions &options, DynVec &v, DynVec &r, SolverStats &stats ) const
{
	stats.reset();
 	if( m_data.nContacts() == 0 ) {
		stats.error = 0 ;
		return 0. ;
	}
	double res = solveContact(options, v, r, stats);

#ifndef SILENCE_ARGUS
	std::cout << "ClothFriction: res = " << res << " \t time=" << stats.time << "s" << std::endl;
#endif

	return res ;
}

double ClothFrictionSolver::solveContact(const SolverOptions &options, DynVec &v, DynVec &r, SolverStats &stats ) const
{
	 if( options.algorithm == SolverOptions::NodalContact ) {
		 return NodalContact( m_data ) .solve( options, v, r, stats );
	 } else if( options.algorithm == SolverOptions::NodalSelfContact) {
		 return NodalSelfContact( m_data ) .solve( options, v, r, stats );
	 } else if (options.algorithm == SolverOptions::TotalContact) {
		return TotalContact(m_data) .solve(options, v, r, stats);
	 }

	 bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > H ;
	 DynVec w ;
	 m_data.makeRelativeVelocityAffineMap( H, w );

	 if( options.algorithm == SolverOptions::AlternatingMinimization ) {
		 return solveAMA( m_data, H, w, options, v, stats ) ;
	 }

	 if( options.algorithm == SolverOptions::ICAGaussSeidel ) {
		 return ICAFrictionSolver(m_data, H, w).solve( options, v, r, stats ) ;
	 }

	 DualFrictionSolver dualSolver( m_data, H, w ) ;
	 return dualSolver.solve( options, v, r, stats ) ;
}


} //argus
