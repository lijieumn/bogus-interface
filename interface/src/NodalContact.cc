
#include "NodalContact.hh"

#include "ClothFrictionData.hh"
#include "SolverOptions.hh"
#include "SolverStats.hh"
#include "RevCoulombLaw.hh"
#include "NodalGaussSeidel.impl.hh"

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>
#include <bogus/Extra/SecondOrder.impl.hpp>

#include <bogus/Interfaces/Cadoux.hpp>
#include <bogus/Interfaces/FrictionProblem.hpp>

#include <bogus/Core/Utils/Timer.hpp>

namespace argus {

double NodalContact::solve(const SolverOptions &options, DynVec& v, SolverStats &stats) const
{
	DynVec r;
	return solve(options, v, r, stats);
}

double NodalContact::solve(const SolverOptions &options, DynVec& v, DynVec& r, SolverStats &stats) const
{
	const int nVertices = m_data.M.rowsOfBlocks() ;
	const int nContacts = m_data.nContacts() ;

	// Convention : set mu=-1 for vertices that are not in contact
	DynVec mu = -DynVec::Ones(nVertices) ;

	bogus::SparseBlockMatrix< Mat3> G ;

	// Fill matrix G with identity matrics
	G.clear();
	G.reserve(nVertices);
	G.setRows(nVertices);
	G.setCols(nVertices);

	for( int i = 0 ; i < nVertices ; ++ i )
	{
		G.insertBack(i,i).setIdentity() ;
	}
	G.finalize() ;

	// Contact basis and external object velocities
	DynVec objVel = DynVec::Zero( m_data.M.rows() ) ;

	for(int i=0; i<nContacts; ++i) {
		int k;
		m_data.coordsA.col(i).maxCoeff(&k);  //Find closest vertex to contact point
		const int vertexId = m_data.indicesA( k, i ) ;

		G.block(vertexId) = m_data.contactBasis.diagonal(i).transpose();

		if( m_data.indicesB(0,i) == -1 ){
			// Compute external object velocity in contact basis
			objVel.segment<3>(vertexId*3) = - G.block(vertexId) * m_data.affineVel.col(i) ;
		}

		mu[vertexId] = options.maxOuterIterations == 0
		        ? m_data.mu[i]                             // Direct GS
		        : 1./(options.tolerance + m_data.mu[i] );  // Reverse Cadoux
	}

	// Pin constraints for duplicated nodes

	typedef bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > PinMat ;
	PinMat P ;
	P.setRows( nVertices ) ;
	P.setCols( nVertices ) ;
	for( unsigned i = 0 ; i < m_data.duplicatedNodes.size() ; ++i )
	{
		const unsigned nDup = m_data.duplicatedNodes[i].size() ;
		for( unsigned j = 0 ; j < nDup ; ++j ) {
			for( unsigned k = 0 ; k < nDup ; ++k ) {
				if( j == k) {
					P.insert( m_data.duplicatedNodes[i][j], m_data.duplicatedNodes[i][j] ) =  Mat3::Identity()/nDup * (nDup - 1.) ;
				} else {
					P.insert( m_data.duplicatedNodes[i][j], m_data.duplicatedNodes[i][k] ) = -Mat3::Identity()/nDup ;
				}
			}
		}
	}
	P.finalize() ;

	typedef bogus::SparseBlockMatrix< Mat3 > ConstraintMat ;
	ConstraintMat C = G * P * G.transpose() ;

	// Prepare system r=Wu+b for Gauss-Seidel solver

	typedef bogus::DualFrictionProblem< 3u >::WType WType ;
	WType W = G * m_data.M * G.transpose() ;
	W.prune(1.e-17) ;


	DynVec b = - G * m_data.f + W * objVel ;

	// Compute initial guess for relative velocities
	DynVec u = G * v - objVel ;

	// Set-up gauss-seidel solver
	NodalGaussSeidel< WType, ConstraintMat > gs( W, C, options.m_inv ) ;
	gs.setTol( options.tolerance );
	gs.setMaxIters( options.maxIterations );
	gs.useInfinityNorm( options.useInfinityNorm );
	gs.setConstraintStepSize( options.nodalConstraintStepSize );
	gs.setConstraintIterations( options.nodalConstraintIterations );
	gs.coloring().update( options.useColoring, W );

	// Solve using Cadoux fixed-point algorithm
	bogus::Signal<unsigned, double> callback ;
	callback.connect( stats, &SolverStats::ack ) ;

	double res = -1 ;

	if( options.maxOuterIterations > 0 ) {
		res = argus::solveCadouxVelWithConstraints< 3 >
		        ( W, b.data(), mu.data(),
		          gs, u.data(),
		          options.maxOuterIterations, &callback ) ;
	} else {
		RevCoulomb3D law( mu.size(), mu.data() ) ;
		law.setLocalTolerance( options.nodalGSLocalTolerance );
		gs.setAutoRegularization( options.nodalGSAutoRegularization );
		gs.callback().connect( callback );
		res = gs.solve( law, b, u ) ;
	}

	// Compute world vel from local vel
	v = G.transpose() * ( u + objVel ) ;
	r = G.transpose() * ( W*u + b );

	return res ;
}

} // argus
