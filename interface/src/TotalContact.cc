#include "TotalContact.hh"

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

double TotalContact::solve(const SolverOptions& options, DynVec &v, SolverStats& stats ) const
{
	DynVec r;
	return solve(options, v, r, stats);
}

void outputMatrix(ClothFrictionData::StiffnessMatrixType M, string name) {
	ofstream out(("make" + name + ".m").c_str());
	out.precision(15);
	out << name << " = [" << std::endl;
	for( int row = 0 ; row < M.rowsOfBlocks() ; ++ row ) {
		double **rowBlock = new double*[3];
		for (int i = 0; i < 3; i++) {
			rowBlock[i] = new double[M.rows()];
			for (int j = 0; j < M.rows(); j++) {
				rowBlock[i][j] = 0;
			}
		}
		for( ClothFrictionData::StiffnessMatrixType::InnerIterator it ( M.innerIterator( row ) ) ; it ; ++it ) {
			for (int i = 0; i < 9; i++) {
				Mat3 b = M.block(it.ptr());
				rowBlock[i/3][it.inner()*3 + i%3] = b(i/3, i%3);
			}
		}
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < M.rows(); j++) {
				out << "\t" << rowBlock[i][j];
			}
			out << ";" << std::endl;
		}
	}
	out << "];" << std::endl;
}

double TotalContact::solve(const SolverOptions &options, DynVec& v, DynVec &r, SolverStats &stats) const
{
	const int nVertices = m_data.M.rowsOfBlocks() ;
	const int nContacts = m_data.nContacts() ;
	std::cout << "in TotalContact solver" << std::endl;
	std::cout << "nVertices:" << nVertices << std::endl;
	std::cout << "nContacts:" << nContacts << std::endl;

	// Convention : set mu=-1 for vertices that are not in contact
	DynVec mu = -DynVec::Ones(nVertices) ;

	bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > G ;
	G.clear();
	G.reserve(nVertices);
	G.setRows(nVertices);
	G.setCols(nVertices);

	// Contact basis and external object velocities
	DynVec objVel = DynVec::Zero( m_data.M.rows() ) ;

	const double eps = 1.e-6 ;
	std::vector< bool > hasContact( nVertices, false ) ;
	double scale = 1;
	for(int i=0; i<nContacts; ++i) {

		vector<int> ci;	// contact indices
		vector<double> cc;	// contact coordinates
		for (int k = 0; k < 3; k++) {
			if( m_data.coordsA(k,i) > eps) {
				ci.push_back(m_data.indicesA(k, i));
				cc.push_back(m_data.coordsA(k, i));
			}
			if(m_data.indicesB(k,i) >= 0 && m_data.coordsB(k,i) > eps ) {
				ci.push_back(m_data.indicesB(k, i));
				cc.push_back(m_data.coordsB(k, i));
			}
		}
		Eigen::VectorXd* rows = NULL;
		if (ci.size() == 4) {
			rows = new Eigen::VectorXd[4];
			for (int r = 0; r < 4; r++) {
				rows[r].resize(4);
			}
			rows[0] << cc[0], cc[1], cc[2], cc[3];
			scale = rows[0].norm();
			rows[0] = rows[0].normalized();
			rows[1] << -cc[1], cc[0], 0, 0;
			rows[1] = rows[1].normalized();
			rows[2] << -cc[2], 0, cc[0], 0;
			rows[3] << -cc[3], 0, 0, cc[0];
			for (int k = 2; k <= 3; k++) {
				rows[k] = rows[k] - rows[k].dot(rows[k - 1])*rows[k - 1];
				rows[k] = rows[k].normalized();
			}
		} else if (ci.size() == 3) {
			rows = new Eigen::VectorXd[3];
			for (int r = 0; r < 3; r++) {
				rows[r].resize(3);
			}
			rows[0] << cc[0], cc[1], cc[2];
			scale = rows[0].norm();
			rows[0] = rows[0].normalized();
			rows[1] << -cc[1], cc[0], 0;
			rows[1] = rows[1].normalized();
			rows[2] << -cc[2], 0, cc[0];
			rows[2] = rows[2] - rows[2].dot(rows[1])*rows[1];
			rows[2] = rows[2].normalized();
		} else if (ci.size() == 2) {
			rows = new Eigen::VectorXd[2];
			for (int r = 0; r < 2; r++) {
				rows[r].resize(2);
			}
			rows[0] << cc[0], cc[1];
			scale = rows[0].norm();
			rows[0] = rows[0].normalized();
			rows[1] << cc[1], -cc[0];
			rows[1] = rows[1].normalized();
		}
		if (ci.size() == 1) {
			hasContact[ci[0]] = true;
			G.insert(ci[0], ci[0]) = m_data.contactBasis.diagonal(i).transpose();
		} else {
			for (unsigned int bi = 0; bi < ci.size(); bi++) {
				hasContact[ci[bi]] = true;
				for (unsigned int bj = 0; bj < ci.size(); bj++) {
					G.insert(ci[bi], ci[bj]) = rows[bi][bj]*m_data.contactBasis.diagonal(i).transpose();
				}
			}
		}
		objVel.segment<3>(ci[0]*3) = - 1/scale * m_data.contactBasis.diagonal(i).transpose() * m_data.affineVel.col(i) ;

		// int k;
		// m_data.coordsA.col(i).maxCoeff(&k);  //Find closest vertex to contact point
		// const int vertexId = m_data.indicesA( k, i ) ;
		// m_data.coordsB.col(i).maxCoeff(&k);  //Find closest vertex to contact point
		// const int vertexJd = m_data.indicesB( k, i );

		// hasContact[vertexId] = true ;

		// // Self-contact case
		// if( vertexJd >= 0 ){
		// 	hasContact[vertexJd] = true ;
		// 	G.insert(vertexId,vertexId) = 1.0 / sqrt(2.0) * m_data.contactBasis.diagonal(i).transpose();
		// 	G.insert(vertexId,vertexJd) = 1.0 / sqrt(2.0) * -m_data.contactBasis.diagonal(i).transpose();
		// 	G.insert(vertexJd,vertexId) = 1.0 / sqrt(2.0) * m_data.contactBasis.diagonal(i).transpose();
		// 	G.insert(vertexJd,vertexJd) = 1.0 / sqrt(2.0) * m_data.contactBasis.diagonal(i).transpose();

		// // Not self-contact case (colliding with external obstacle)
		// } else {
		// 	G.insert(vertexId,vertexId) = m_data.contactBasis.diagonal(i).transpose();
		// }
		// objVel.segment<3>(vertexId*3) = - m_data.contactBasis.diagonal(i).transpose() * m_data.affineVel.col(i) ;
		mu[ci[0]] = options.maxOuterIterations == 0
		        ? m_data.mu[i]                             // Direct GS
		        : 1./(options.tolerance + m_data.mu[i] );  // Reverse Cadoux
	}
	for( int i = 0 ; i < nVertices ; ++ i )
	{
		if(!hasContact[i])
			G.insertBack(i,i).setIdentity() ;
	}

	G.finalize();
	std::cout << "number of blocks:" << G.nBlocks() << std::endl;

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

	typedef bogus::SparseBlockMatrix< Mat3 > ConstraintMat ;
	ConstraintMat C = G * P * G.transpose() ;

	// Prepare system r=Wu+b for Gauss-Seidel solver

	typedef bogus::DualFrictionProblem< 3u >::WType WType ;
	WType W = G * m_data.M * G.transpose() ;
	W.prune(1.e-17) ;
	// outputMatrix(m_data.M, "M");
	// outputMatrix(G, "G");

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
	r = G.transpose() * (W*u + b);

	return res ;
}

} // argus
