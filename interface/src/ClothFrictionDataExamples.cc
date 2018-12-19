#include "ClothFrictionData.hh"

#include <Eigen/Geometry>

#include <bogus/Core/Block.impl.hpp>

namespace argus {

static void makeNormalMatrix( const Eigen::Vector3d &n, Eigen::Matrix3d& mat)
{
	mat.col(0) = n ;
	if(n[0] > 0) {
		mat.col(1) = Eigen::Vector3d(-n[1], n[0], 0).normalized() ;
	} else {
		mat.col(1) = Eigen::Vector3d(0, -n[2], n[1]).normalized() ;
	}
	mat.col(2) = mat.col(0).cross( mat.col(1) ) ;
}



bool ClothFrictionData::loadSimpleTest(){

	// We create a simple test where a triangle collides with an external object

	// Stiffness matrix //
	M.clear();
	M.reserve(3);
	M.setRows(3);
	M.setCols(3);
	M.insertBack(0,0).setIdentity();   // vert 0
	M.insertBack(1,1).setIdentity();   // vert 1
	M.insertBack(2,2).setIdentity();   // vert 2
	M.finalize();

	// Forces on each of the three triangle vertices
	f.resize(9);
	f.segment(0,3) = Eigen::Vector3d(9.8,0,0);
	f.segment(3,3) = Eigen::Vector3d(9.8,0,0);
	f.segment(6,3) = Eigen::Vector3d(9.8,0,0);

	// Friction coefficients for each contact //
	mu.resize(2);
	mu[0] = 1;

	// Contact Basis //
	contactBasis.clear();
	contactBasis.reserve(1);
	contactBasis.setRows(1);
	contactBasis.setCols(1);
	contactBasis.insertBack(0,0).setIdentity();  // because vertex 0 is coming in from -x direction, normal is +x and contact matrix is identity
	contactBasis.finalize() ;

	//std::cout << "Contact basis: \n";
	//std::cout << "C00: \n" << contactBasis.block( contactBasis.blockPtr(0,0) ) << std::endl;
	//std::cout << "C11: \n" << contactBasis.block( contactBasis.blockPtr(1,1) ) << std::endl;
	//std::cout << "C22: \n" << contactBasis.block( contactBasis.blockPtr(2,2) ) << std::endl;


	// Three indices and corresponding coords for the triangle
	indicesA.resize(3,1);
	coordsA.resize(3,1);
	indicesB.resize(3,1);
	coordsB.resize(3,1);
	affineVel.resize(3,1);

	indicesA.col(0) = Eigen::Vector3i(0,1,2);
	coordsA.col(0) = Eigen::Vector3d(1,0,0);  // collision occurs exactly at vertex 0
	indicesB.col(0) = Eigen::Vector3i(-1,-1,-1);
	coordsB.col(0).setZero();
	affineVel.col(0) = Eigen::Vector3d(0,-1,0);  // assume external object at rest

	/*indicesA.col(1) = Eigen::Vector3i(0,1,2);
	coordsA.col(1) = Eigen::Vector3d(0,1,0);  // collision occurs exactly at vertex 0
	indicesB.col(1) = Eigen::Vector3i(-1,-1,-1);
	coordsB.col(1).setZero();
	affineVel.col(1) = Eigen::Vector3d(0,0,0);  // assume external object at rest  */

	return true;

}

bool ClothFrictionData::loadSlidingTest(){

	// We create a simple test where a triangle collides with an external object

	// Stiffness matrix //
	M.clear();
	M.reserve(3);
	M.setRows(3);
	M.setCols(3);
	M.insertBack(0,0).setIdentity();   // vert 0
	M.insertBack(1,1).setIdentity();   // vert 0
	M.insertBack(2,2).setIdentity();   // vert 0
	M.finalize();

	// Forces on each of the three triangle vertices
	f.resize(9);
	f.segment(0,3) = Eigen::Vector3d(-4,0,-9.8);
	f.segment(3,3) = Eigen::Vector3d(-4,0,-9.8);
	f.segment(6,3) = Eigen::Vector3d(-4,0,-9.8);

	// Friction coefficients for each of the three triangle vertices //
	mu.resize(3);
	mu[0] = 1;
	mu[1] = 1;
	mu[2] = 1;

	// Contact Basis //
	contactBasis.clear();
	contactBasis.reserve(3);
	contactBasis.setRows(3);
	contactBasis.setCols(3);
	// Slightly tilted normals
	makeNormalMatrix( Eigen::Vector3d(1,-.1,-.1).normalized(), contactBasis.insertBack(0,0));
	makeNormalMatrix( Eigen::Vector3d(1,.2,0).normalized(), contactBasis.insertBack(1,1));
	makeNormalMatrix( Eigen::Vector3d(1,0,.2).normalized(), contactBasis.insertBack(2,2));
	contactBasis.finalize() ;

	// Three indices and corresponding coords for the triangle
	indicesA.resize(3,3);
	coordsA.resize(3,3);
	indicesB.resize(3,3);
	coordsB.resize(3,3);
	affineVel.resize(3,3);

	affineVel.col(0) = Eigen::Vector3d(-1,0,0);
	affineVel.col(1) = Eigen::Vector3d(-1,0,0);
	affineVel.col(2) = Eigen::Vector3d(-1,0,0);

	indicesA.col(0) = Eigen::Vector3i(0,1,2);
	indicesA.col(1) = Eigen::Vector3i(0,1,2);
	indicesA.col(2) = Eigen::Vector3i(0,1,2);


	coordsA.col(0) = Eigen::Vector3d(1,0,0);  // collision occurs exactly at vertex 0
	coordsA.col(1) = Eigen::Vector3d(1,0,0);  // collision occurs exactly at vertex 0
	coordsA.col(2) = Eigen::Vector3d(1,0,0);  // collision occurs exactly at vertex 0

	indicesB.col(0) = Eigen::Vector3i(-1,-1,-1);
	indicesB.col(1) = Eigen::Vector3i(-1,-1,-1);
	indicesB.col(2) = Eigen::Vector3i(-1,-1,-1);

	coordsB.setZero();

	return true;

}

bool ClothFrictionData::loadPinTest(){
	ClothFrictionData orig ;
	orig.loadSimpleTest() ;

	const int I = 0 ; // index of vertex to duplicate
	const int n = orig.M.rowsOfBlocks() ;

	StiffnessMatrixType S ; //scaling-duplication matrix
	S.setRows( n+1 ) ;
	S.setCols( n ) ;

	for( int i = 0, j = 0 ; i < orig.M.rowsOfBlocks() ; ++i, ++j ) {
		if( i == I ) {
			S.insertBack( j, i ) = Mat3::Identity() * .5 ;
			++j ;
			S.insertBack( j, i ) = Mat3::Identity() * .5 ;
		} else {
			S.insertBack( j, i ).setIdentity() ;
		}
	}
	S.finalize() ;

	M = S * orig.M * S.transpose() ;

	// Unsingularize
	M.block(I,I) = .5 * orig.M.block(I,I) ;
	M.block(I,I+1).setZero() ;
	M.block(I+1,I).setZero() ;
	M.block(I+1,I+1) = .5 * orig.M.block(I,I) ;

	f = S * orig.f ;

	mu = orig.mu ;

	contactBasis = orig.contactBasis ;

	indicesA = orig.indicesA ;
	coordsA = orig.coordsA ;
	indicesB = orig.indicesB ;
	coordsB = orig.coordsB ;

	for( int i = 0 ; i < nContacts() ; ++i ) {
		for( int j = 0 ; j < 3 ; ++j ) {
			if( indicesA(j,i) > I ) {
				indicesA(j,i) += 1 ;
			}
		}
	}

	affineVel = orig.affineVel ;

	std::vector<int> dup ;
	dup.push_back(I) ;
	dup.push_back(I+1) ;
	duplicatedNodes.push_back(dup);

	return true ;
}

bool ClothFrictionData::loadTwoContactsTest(){


	// We create a simple test where a triangle collides with two external objects

	// Stiffness matrix //
	M.clear();
	M.reserve(3);
	M.setRows(3);
	M.setCols(3);
	M.insertBack(0,0).setIdentity();   // vert 0
	M.insertBack(1,1).setIdentity();   // vert 1
	M.insertBack(2,2).setIdentity();   // vert 2
	M.finalize();

	// Forces on each of the three triangle vertices
	f.resize(9);
	f.segment(0,3) = Eigen::Vector3d(0,0,-9.8);
	f.segment(3,3) = Eigen::Vector3d(0,0,-9.8);
	f.segment(6,3) = Eigen::Vector3d(0,0,-9.8);

	// Friction coefficients for each of the three triangle vertices //
	mu.resize(2);
	mu[0] = 0.1;
	mu[1] = 0.1;

	// Contact Basis //
	int numContact = 2;
	contactBasis.clear();
	contactBasis.reserve(numContact);
	contactBasis.setRows(numContact);
	contactBasis.setCols(numContact);
	contactBasis.insertBack(0,0).setIdentity();
	contactBasis.insertBack(1,1) << 0, 1, 0, 0, 0, 1, 1, 0, 0;
	contactBasis.finalize() ;

	// Three indices and corresponding coords for the triangle
	indicesA.resize(3,numContact);
	coordsA.resize(3,numContact);
	indicesB.resize(3,numContact);
	coordsB.resize(3,numContact);
	affineVel.resize(3,numContact);

	indicesA.col(0) = Eigen::Vector3i(0,1,2);
	coordsA.col(0) = Eigen::Vector3d(0,1,0);  // collision occurs exactly at vertex 1

	indicesB.col(0) = Eigen::Vector3i(-1,-1,-1);
	coordsB.col(0).setZero();

	indicesA.col(1) = indicesA.col(0) ;
	indicesB.col(1) = indicesB.col(0) ;
	coordsA.col(1) = coordsA.col(0) ;
	coordsB.col(1) = coordsB.col(0) ;

	affineVel.col(0) = Eigen::Vector3d(-1,0.3, 0);  // assume external object at rest
	affineVel.col(1) = Eigen::Vector3d(-1,0,-0.2);  // assume external object at rest

	return true;

}

bool ClothFrictionData::loadTwoLayersTest(){


	// We create a simple test where two triangles are on top of  an external object

	// Stiffness matrix //
	M.clear();
	M.reserve(4);
	M.setRows(4);
	M.setCols(4);
	M.insertBack(0,0).setIdentity();   // vert 0
	M.insertBack(1,1).setIdentity();   // vert 1
	M.insertBack(2,2).setIdentity();   // vert 2
	M.insertBack(3,3).setIdentity();   // vert 3
	M.finalize();

	// Forces on each of the three triangle vertices
	f.resize(12);
	f.segment(0,3) = Eigen::Vector3d(0,0,-9.8);
	f.segment(3,3) = Eigen::Vector3d(0,0,-9.8);
	f.segment(6,3) = Eigen::Vector3d(0,0,-9.8);
	f.segment(9,3) = Eigen::Vector3d(0,0,-9.8);

	// Friction coefficients for each of the three triangle vertices //
	mu.resize(3);
	mu[0] = 0.5;
	mu[1] = 0.2;
	mu[2] = 0.;

	// Contact Basis //
	contactBasis.clear();
	contactBasis.reserve(3);
	contactBasis.setRows(3);
	contactBasis.setCols(3);
	contactBasis.insertBack(0,0) << 0, 1, 0, 0, 0, 1, 1, 0, 0;
	contactBasis.insertBack(1,1) << 0, 0, -1, 0, -1, 0, -1, 0, 0;
	contactBasis.insertBack(2,2).setIdentity() ;
	contactBasis.finalize() ;

	// Three indices and corresponding coords for the triangle
	indicesA.resize(3,3);
	coordsA.resize(3,3);
	indicesB.resize(3,3);
	coordsB.resize(3,3);
	affineVel.resize(3,3);

	indicesA.col(0) = Eigen::Vector3i(0,1,2);
	coordsA.col(0) = Eigen::Vector3d(0,1,0);  // collision occurs exactly at vertex 1
	indicesB.col(0) = Eigen::Vector3i(-1,-1,-1);
	coordsB.col(0).setZero();

	indicesA.col(1) = indicesA.col(0) ;
	coordsA.col(1) = coordsA.col(0) ;
	indicesB.col(1) = Eigen::Vector3i(2,3,0);
	coordsB.col(1) = Eigen::Vector3d(0,1,0) ;
	//indicesB.col(1) = Eigen::Vector3i(-1,-1,-1);
	//coordsB.col(1).setZero();

	indicesA.col(2) = indicesA.col(0) ;
	coordsA.col(2) = coordsA.col(0) ;
	indicesB.col(2) = indicesB.col(0) ;
	coordsB.col(2) = coordsB.col(0) ;

	affineVel.col(0) = Eigen::Vector3d(-1,0.3, 0);
	affineVel.col(1).setZero();
	affineVel.col(2) = Eigen::Vector3d(-1.2,0.3, 0);

	return true;

}


bool ClothFrictionData::loadSelfIntersectTest(){

	// We create a simple test where two verts from the same mesh collide with each other

	// Stiffness matrix //
	M.clear();
	M.reserve(2);
	M.setRows(2);
	M.setCols(2);
	M.insertBack(0,0).setIdentity();   // vert 0 (tri A)
	M.insertBack(1,1).setIdentity();   // vert 0 (tri B)
	M.finalize();
	M.block( M.blockPtr(0,0) )(0,0) *= 2;
	M.block( M.blockPtr(0,0) )(1,1) *= 2;
	M.block( M.blockPtr(0,0) )(2,2) *= 2;

	std::cout << "Stiffness matrix: \n";
	std::cout << "M00: \n" << M.block( M.blockPtr(0,0) ) << std::endl;
	std::cout << "M01: \n" << M.block( M.blockPtr(0,1) ) << std::endl;
	std::cout << "M10: \n" << M.block( M.blockPtr(0,1) ) << std::endl;
	std::cout << "M11: \n" << M.block( M.blockPtr(1,1) ) << std::endl;

	// Forces on each of the verts
	f.resize(6);
	f.segment(0,3) = Eigen::Vector3d(-10,0,0);  // force on vert 0
	f.segment(3,3) = Eigen::Vector3d(0,0,0); // force on vert 3

	// Friction coefficient //
	mu.resize(1);
	mu[0] = 1.0;

	// Contact Basis //
	// Vert 0 is coming from the -x direction, normal is +x and contact matrix is identity
	// Vert 1 is coming from the +x direction, normal is -x and contact matrix is diag(-1,-1,1)
	contactBasis.clear();
	contactBasis.reserve(2);
	contactBasis.setRows(2);
	contactBasis.setCols(2);
	contactBasis.insertBack(0,0).setIdentity();
	contactBasis.insertBack(0,1).setZero();
	contactBasis.insertBack(1,0).setZero();
	contactBasis.insertBack(1,1).setIdentity();
	contactBasis.finalize();
	contactBasis.block( contactBasis.blockPtr(1,1) )(0,0) *= 1;
	contactBasis.block( contactBasis.blockPtr(1,1) )(1,1) *= 1;
	contactBasis.block( contactBasis.blockPtr(1,1) )(2,2) *= 1;

	std::cout << "Contact basis: \n";
	std::cout << "C00: \n" << contactBasis.block( contactBasis.blockPtr(0,0) ) << std::endl;
	std::cout << "C01: \n" << contactBasis.block( contactBasis.blockPtr(0,1) ) << std::endl;
	std::cout << "C10: \n" << contactBasis.block( contactBasis.blockPtr(1,0) ) << std::endl;
	std::cout << "C11: \n" << contactBasis.block( contactBasis.blockPtr(1,1) ) << std::endl;

	// Three indices and corresponding coords for the triangle
	indicesA.resize(3,1);
	coordsA.resize(3,1);
	indicesB.resize(3,1);
	coordsB.resize(3,1);
	affineVel.resize(3,1);

	indicesA.col(0) = Eigen::Vector3i(0,2,3);
	coordsA.col(0) = Eigen::Vector3d(1,0,0);  // collision occurs exactly at vertex 0
	indicesB.col(0) = Eigen::Vector3i(1,4,5);
	coordsB.col(0) = Eigen::Vector3d(1,0,0);
	affineVel.col(0) = Eigen::Vector3d(0,0,0);  // assume external object at rest

	return true;


}


} // argus
