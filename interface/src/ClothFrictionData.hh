#ifndef ARGUS_CLOTH_FRICTION_DATA_HH
#define ARGUS_CLOTH_FRICTION_DATA_HH

#include "MatrixTypes.hh"

#include <bogus/Core/Block.hpp>

namespace argus {

//! Representation of a cloth frictional contact problem, with I/O capabilities
struct ClothFrictionData {

	typedef bogus::SparseBlockMatrix< Mat3 > StiffnessMatrixType ;
	typedef bogus::SparseBlockMatrix< Mat3 > ContactBasisMatrix ;

	// Dynamics

	StiffnessMatrixType M ;  //!< Stiffness matrix
	DynVec				f ;  //!< Forces, such that the unconstrained dynamics are the solutuon of  Mv = f

	StiffnessMatrixType S ; //scaling-duplication matrix
	StiffnessMatrixType pinvS ; //psudo inverse of scaling-duplication matrix

	// Contacts

	DynVec  mu  ; 		    //!< Friction coefficient at each contact
	//! Bloc diagonal matrix of rotations from local contact basis to world coordinates
	/*! That is, each 3x3 diagonal block of contactBasis should be (n|t1|t2), where n is the contac normal
	 *  and t1 and t2 are two unite orthogonal tangent vectors
	 */
	ContactBasisMatrix	contactBasis ;

	DynMat3i  indicesA; //!< Indices of the vertices of the first contacting triangle, one column per contact
	DynMat3    coordsA; //!< Barycentric coordinates of the first contacting triangle, one column per contact

	DynMat3i  indicesB; //!< Indices of the vertices of the second contacting triangle (or -1 for external object)
	DynMat3    coordsB; //!< Barycentric coordinates of the second contacting triangle

	//! Affine part of the relative velocity computation.
	/*! For cloth/mesh contacts, should be the negative of the object's velocity
	 *  For continous-time collision detection, one may add the offset (xA-xB) /dt
	 *  to account for the displacement of the two objects before the collision.
B	 *  where xZ is the position of the Z object, and n the contact normal (contactBasis.col(3*i))
	 */
	DynMat3  affineVel;

	std::vector< std::vector< int > > duplicatedNodes ; //! List of (sets  of node indices that should be pinned together) (NodalContact only)

	int nContacts() const { return indicesA.cols() ; } //!< Number of contacts
	int nVertices() const { return M.rowsOfBlocks() ; } //!< Number of vertices (including duplicated)
	int nDuplicatedVertices() const {
		int ndv = 0 ;
		for( unsigned i = 0 ; i < duplicatedNodes.size(); ++i )
			ndv += duplicatedNodes[i].size() - 1 ;
		return ndv ;
	}

	// Adhesion
	void applyAdhesion( const DynVec& adh ) ;

	void alignTangentsWithVelocity(const DynVec &v) ;

	void makeRelativeVelocityAffineMap (
	    bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > & H,
	    DynVec &w, bool snapToVertices = false ) const ;

	// Duplication

	//! Finds nodes with several contacts and duplicates them
	bool findAndDuplicate() ;
	//! Copies velocity of original nodes to that of duplicated nodes
	void scatterVelocities( const DynVec& orig, DynVec& duplicated ) const ;
	//! Merges velocities of duplicated nodes
	void  gatherVelocities( const DynVec& duplicated, DynVec& merged ) const ;

	// I/O

	template <typename Archive>
	void serialize(Archive &ar, const unsigned int ) ;

	bool dump( const std::string &file ) const ;
	bool load( const std::string &file ) ;
	bool loadSimpleTest();
	bool loadPinTest();
	bool loadTwoContactsTest() ;
	bool loadTwoLayersTest() ;
	bool loadSlidingTest() ;
	bool loadSelfIntersectTest();
};


} //argus

#endif
