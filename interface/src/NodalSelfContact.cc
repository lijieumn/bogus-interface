#include "NodalSelfContact.hh"

#include "ClothFrictionData.hh"
#include "SolverOptions.hh"
#include "SolverStats.hh"
#include "RevCoulombLaw.hh"
#include "NodalGaussSeidel.impl.hh"
#include <fstream>
using namespace std;

#include <bogus/Core/Block.impl.hpp>
#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/BlockSolvers/GaussSeidel.impl.hpp>
#include <bogus/Extra/SecondOrder.impl.hpp>

#include <bogus/Interfaces/Cadoux.hpp>
#include <bogus/Interfaces/FrictionProblem.hpp>

#include <bogus/Core/Utils/Timer.hpp>

namespace argus {

typedef bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > SBM;
void printMatrix(SBM M, string name) {
	int n = M.rowsOfBlocks();
	string file = "make" + name + ".m";
	ofstream out(file);
	double row[3][n*3];
	out << name << " = [" << std::endl;
	for (int i = 0; i < n; i++) {
		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3*n; jj++) {
				row[ii][jj] = 0;
			}
		}
		for (SBM::InnerIterator it (M.innerIterator(i)); it ; ++ it) {
			int col = it.inner();
			Mat3 b = M.block(it.ptr());
			for (int k = 0; k < 9; k++) {
				row[k/3][col*3 + k%3] = b(k/3, k%3);
			}
		}
		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3*n; jj++) {
				out << row[ii][jj] << "\t";
			}
			out << std::endl;
		}
	}
	out << "]";
}

struct NodalOverrideFunctor : public EvalFunctor
{

	typedef bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > GType ;
	typedef bogus::DualFrictionProblem< 3u >::WType WType ;
	typedef bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > PinMat ;

	NodalOverrideFunctor(
	        std::unique_ptr<EvalFunctor> && orig ,
	        const GType& G,
	        const PinMat& P,
	        const DynVec &objVel,
	        const DynVec &b,
	        std::vector<int>& contactIds,
	        DynVec& u,
	        DynVec& v,
	        DynVec& rLoc,
	        DynVec& rGlob
	        )
	    : m_orig(std::move(orig))
	    , m_G(G)
	    , m_P(P)
	    , m_objVel(objVel)
	    , m_b(b)
	    , m_contactIds(contactIds)
	    , m_u(u)
	    , m_v(v)
	    , m_rLoc(rLoc)
	    , m_rGlob(rGlob)
	{

	}

	double operator() () const override
	{
		m_v = m_G.transpose() * ( m_u + m_objVel ) ;

		for( unsigned vi = 0 ; vi < m_contactIds.size() ; ++vi )
		{
			if( -1 != m_contactIds[vi]) {
				m_rGlob.segment<3>(3*m_contactIds[vi]) = m_rLoc.segment<3>(3*vi);
			}
		}

		return (*m_orig)();
	}

	std::unique_ptr<EvalFunctor> m_orig ;

	const GType&  m_G ;
	const PinMat& m_P ;

	const DynVec &m_objVel;
	const DynVec &m_b;
	const std::vector<int>& m_contactIds;

	DynVec& m_u ;
	DynVec& m_v;
	DynVec& m_rLoc;
	DynVec& m_rGlob;
};


double NodalSelfContact::solve(const SolverOptions& options, DynVec &v, SolverStats& stats ) const
{
	DynVec r;
	return solve(options, v, r, stats);
}

double NodalSelfContact::solve(const SolverOptions &options, DynVec& v, DynVec &r, SolverStats &stats) const
{
	const int nVertices = m_data.M.rowsOfBlocks() ;
	const int nContacts = m_data.nContacts() ;

	// Convention : set mu=-1 for vertices that are not in contact
	DynVec mu = -DynVec::Ones(nVertices) ;


	typedef bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > GType ;
	GType G;
	G.clear();
	// G.reserve(nVertices);
	G.setRows(nVertices);
	G.setCols(nVertices);

	// Contact basis and external object velocities
	DynVec objVel = DynVec::Zero( m_data.M.rows() ) ;

	std::vector< bool > hasContact( nVertices, 0 ) ;
	std::vector< int  > contactIds( nVertices, -1 ) ;
	for(int i=0; i<nContacts; ++i) {

		int k;
		m_data.coordsA.col(i).maxCoeff(&k);  //Find closest vertex to contact point
		const int vertexId = m_data.indicesA( k, i ) ;
		m_data.coordsB.col(i).maxCoeff(&k);  //Find closest vertex to contact point
		const int vertexJd = m_data.indicesB( k, i );

		contactIds[vertexId] = true ;

		const Mat3& Ei = m_data.contactBasis.diagonal(i).transpose() ;

		// Self-contact case
		if( vertexJd >= 0 ){
			hasContact[vertexJd] = true;

			G.insert(vertexId,vertexId) = 1.0 / sqrt(2.0) * Ei ;
			G.insert(vertexId,vertexJd) = 1.0 / sqrt(2.0) * -Ei;
			G.insert(vertexJd,vertexId) = 1.0 / sqrt(2.0) * Ei;
			G.insert(vertexJd,vertexJd) = 1.0 / sqrt(2.0) * Ei;

			objVel.segment<3>(vertexId*3) = - 1.0 / sqrt(2.0) * Ei * m_data.affineVel.col(i) ;
		// Not self-contact case (colliding with external obstacle)
		} else {
			G.insert(vertexId,vertexId) = Ei;
			objVel.segment<3>(vertexId*3) = - Ei * m_data.affineVel.col(i) ;
		}
		contactIds[vertexId] = i;
		hasContact[vertexId] = true;

		mu[vertexId] = options.maxOuterIterations == 0
		        ? m_data.mu[i]                             // Direct GS
		        : 1./(options.tolerance + m_data.mu[i] );  // Reverse Cadoux

	}
	for( int i = 0 ; i < nVertices ; ++ i )
	{
		if( !hasContact[i])
			G.insert(i,i).setIdentity() ;
	}

	G.finalize();

	/*std::cout << "Printing out G matrix" << std::endl;
	for(int i = 0; i < nVertices; ++i){
		for(int j = 0; j < nVertices; ++j){
			if( fabs(G.block( G.blockPtr(i,j))(0,0) ) < 1.e-10
			&&  fabs(G.block( G.blockPtr(i,j))(0,1) ) < 1.e-10
			&&  fabs(G.block( G.blockPtr(i,j))(0,2) ) < 1.e-10
			&&  fabs(G.block( G.blockPtr(i,j))(1,0) ) < 1.e-10
			&&  fabs(G.block( G.blockPtr(i,j))(1,1) ) < 1.e-10
			&&  fabs(G.block( G.blockPtr(i,j))(1,2) ) < 1.e-10
			&&  fabs(G.block( G.blockPtr(i,j))(2,0) ) < 1.e-10
			&&  fabs(G.block( G.blockPtr(i,j))(2,1) ) < 1.e-10
			&&  fabs(G.block( G.blockPtr(i,j))(2,2) ) < 1.e-10 ){
				std::cout << "(" << i << "," << j << "): " << "ZERO" << std::endl;
			} else {
				std::cout << "(" << i << "," << j << "): \n" << G.block( G.blockPtr(i,j) ) << std::endl;
			}
		}
	}*/

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

	typedef bogus::SparseBlockMatrix< Mat3, bogus::SYMMETRIC > ConstraintMat ;
//	typedef bogus::Product< bogus::Product<GType, PinMat>, bogus::Transpose<GType> > ConstraintMat ;
	ConstraintMat C = G * P * G.transpose() ;
	C.cacheTranspose();

	// Prepare system r=Wu+b for Gauss-Seidel solver

	typedef bogus::DualFrictionProblem< 3u >::WType WType ;
	WType W = G * m_data.M * G.transpose() ;
	W.prune(1.e-17) ;
	W.cacheTranspose();
#ifndef SILENCE_ARGUS
	std::cerr << "Num W blocks " << W.nBlocks() << std::endl ;
#endif

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
	gs.setSkipTol( options.tolerance );
//	gs.setSkipIters( 5 );

	gs.setAScaling( options.AScale );
	gs.setBScaling( options.bScale );
	gs.setEvalScaling( options.evalScale );

	gs.coloring().update( options.useColoring, W );


	// Solve using Cadoux fixed-point algorithm
	bogus::Signal<unsigned, double> callback ;
	callback.connect( stats, &SolverStats::ack ) ;

	double res = -1 ;

	DynVec rr (u.rows()) ;
//	r.setZero( 3*nContacts );

	if(stats.functor)
	{
		stats.functor.reset( new NodalOverrideFunctor(
		                         std::move(stats.functor),
		                         G, P, objVel, b, contactIds, u, v, rr, r)) ;
	}

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
		res = gs.solve( law, b, u, rr ) ;
	}

	// std::cout << "in solver:" << std::endl;
	// for( unsigned i = 0 ; i < m_data.duplicatedNodes.size() ; ++i )
	// {
	// 	std::cout << "duplicated:" << std::endl;
	// 	const unsigned nDup = m_data.duplicatedNodes[i].size() ;
	// 	for( unsigned j = 0 ; j < nDup ; ++j ) {
	// 		std::cout << "v:" << v.segment(m_data.duplicatedNodes[i][j]*3, 3) << std::endl;
	// 	}
	// }

	// Compute world vel from local vel
	v = G.transpose() * ( u + objVel ) ;
	r = DynVec( G.transpose() * (W*u + b) );

	return res ;
}

} // argus
