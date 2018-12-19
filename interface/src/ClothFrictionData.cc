#include "ClothFrictionData.hh"

#include <Eigen/Geometry>

#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/Block.impl.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include <boost/version.hpp>

#include <fstream>

BOOST_CLASS_VERSION(argus::ClothFrictionData, 1)

namespace argus {

void ClothFrictionData::scatterVelocities( const DynVec& orig, DynVec& duplicated ) const
{

	const int dVertices = nVertices() ;
	const int oVertices = dVertices - nDuplicatedVertices() ;

//	assert( orig.rows() = 3*oVertices ) ;
	duplicated.resize( 3*dVertices ) ;

	std::vector<bool> taken( dVertices+1, false ) ;

	std::size_t dn_index = 0 ;
	std::size_t cur_dest = 0 ;

	for( int i = 0 ; i < oVertices; ++i ) {

		duplicated.segment<3>(3*cur_dest) = orig.segment<3>(3*i) ;

		if( duplicatedNodes.size() > dn_index
		        && duplicatedNodes[dn_index][0] == (int)cur_dest ) {
			// Duplicated node
			for( unsigned j = 1 ; j < duplicatedNodes[dn_index].size() ; ++j ) {
				const int dup_idx = duplicatedNodes[dn_index][j] ;
				taken[dup_idx] = true ;
				duplicated.segment<3>(3*dup_idx) = orig.segment<3>(3*i) ;
			}
			++ dn_index ;
		}
		// Find next free spot
		while( taken[++cur_dest] ) ;

	}
}


void ClothFrictionData:: gatherVelocities( const DynVec& duplicated, DynVec& merged ) const
{

	const int dVertices = nVertices() ;
	const int oVertices = dVertices - nDuplicatedVertices() ;


	std::vector<bool> taken( dVertices, false ) ;

//	assert( duplicated.rows() = 3*dVertices ) ;
	merged.resize(3*oVertices) ;

	std::size_t dn_index = 0 ;
	std::size_t cur_dest = 0 ;

	for( int i = 0 ; i < dVertices ; ++i ) {
		if( taken[i]) continue ;

		if( duplicatedNodes.size() > dn_index
		        && duplicatedNodes[dn_index][0] == i ) {
			merged.segment<3>(3*cur_dest).setZero() ;
			for( unsigned j = 0 ; j < duplicatedNodes[dn_index].size() ; ++j ) {
				const int dup_idx = duplicatedNodes[dn_index][j] ;
				taken[dup_idx] = true ;
				merged.segment<3>(3*cur_dest) += duplicated.segment<3>(3*dup_idx) ;
			}
			merged.segment<3>(3*cur_dest) /= duplicatedNodes[dn_index].size() ;
			++ dn_index;
		} else {
			merged.segment<3>(3*cur_dest) = duplicated.segment<3>(3*i) ;
		}

		++ cur_dest ;
	}
}

bool ClothFrictionData::findAndDuplicate()
{
	const double eps = 1.e-6 ;

	duplicatedNodes.clear() ;

	const int nVertices = M.rowsOfBlocks() ;
	std::vector<int> pernode(nVertices, 0) ;

	//Count contacts per node
	for( int i = 0 ; i < nContacts() ; ++i ) {

		for( int k = 0 ; k < 3 ; ++k ) {
			if( coordsA(k,i) > eps)
				++pernode[indicesA(k,i)] ;
			if(indicesB(k,i) >= 0 && coordsB(k,i) > eps )
				++pernode[indicesB(k,i)] ;
		}

	}

	// Mark duplicated
	std::vector<int> oindices(nVertices) ;
	std::vector<int> dindices(nVertices, -1) ;

	std::size_t cur_index = 0 ;
	for( int i = 0 ; i < nVertices ; ++i ) {
		oindices[i] = cur_index ;
		if(pernode[i] > 1) {

			duplicatedNodes.push_back(std::vector<int>()) ;
			duplicatedNodes.back().push_back( cur_index ) ;
			for( int j = 1 ; j < pernode[i] ; ++ j )
				duplicatedNodes.back().push_back( ++cur_index ) ;

			dindices[i] = duplicatedNodes.size()-1 ;
		}
		++cur_index ;
	}

	const int dVertices = cur_index;

	// Expand mass matrix and force vector
	// StiffnessMatrixType S ; //scaling-duplication matrix
	S.setRows( dVertices ) ;
	S.setCols( nVertices ) ;
	pinvS.setRows( nVertices );
	pinvS.setCols( dVertices );

	for( int i = 0 ; i < nVertices ; ++i ) {
		if( dindices[i] != -1 ) {
			const double s = 1. / duplicatedNodes[dindices[i]].size() ;
			for( unsigned j = 0 ; j < duplicatedNodes[dindices[i]].size() ; ++ j) {
				S.insertBack( duplicatedNodes[dindices[i]][j], i ) = Mat3::Identity() * s ;
				pinvS.insertBack( i, duplicatedNodes[dindices[i]][j]) = Mat3::Identity() * s;
			}
		} else {
			S.insertBack( oindices[i], i ).setIdentity() ;
			pinvS.insertBack(i, oindices[i]).setIdentity() ;
		}
	}
	S.finalize() ;
	pinvS.finalize();

	{
		StiffnessMatrixType oldM = M ;
		M = S * oldM * S.transpose() ;
		DynVec oldf = f ;
		f = S * oldf ;
	}

	// Unsingularize mass matrix
	for( unsigned i = 0 ; i < duplicatedNodes.size() ; ++i ) {
		for( unsigned j = 0 ; j < duplicatedNodes[i].size() ; ++j ) {
			for( unsigned k = 0 ; k < duplicatedNodes[i].size() ; ++k ) {
				if( j == k )
					M.block(duplicatedNodes[i][j], duplicatedNodes[i][j]) *= duplicatedNodes[i].size() ;
				else
					M.block(duplicatedNodes[i][j], duplicatedNodes[i][k]).setZero() ;
			}
		}
	}

	// Update contact indices
	pernode.assign(nVertices, 0) ;
	for( int i = 0 ; i < nContacts() ; ++i ) {

		for( int k = 0 ; k < 3 ; ++k ) {
			if( coordsA(k,i) > eps ) {
				const int idx = indicesA(k,i) ;
				if( ++pernode[idx] > 1 )
				{
					indicesA(k,i) = duplicatedNodes[dindices[idx]][pernode[idx]-1]  ;
				} else {
					indicesA(k,i) = oindices[idx] ;
				}
			}

			if(indicesB(k,i) >= 0 && coordsB(k,i) > eps ) {
				const int idx = indicesB(k,i) ;
				if( ++pernode[idx] > 1 ) {
					indicesB(k,i) = duplicatedNodes[dindices[idx]][pernode[idx]-1]  ;
				} else {
					indicesB(k,i) = oindices[idx] ;
				}
			}
		}

	}

	return true ;
}


template <class Archive>
void ClothFrictionData::serialize(Archive &ar, const unsigned int version )
{
	ar & M ;
	ar & f ;

	ar & mu ;
	ar & contactBasis ;

	ar & indicesA;
	ar &  coordsA;

	ar & indicesB;
	ar &  coordsB;

	ar & affineVel;

	if(version > 0) {
		ar & duplicatedNodes ;
	}
}

bool ClothFrictionData::load(const std::string &file)
{
	std::ifstream ifs( file.c_str() );

	std::cout << "Boost version: "
	      << BOOST_VERSION / 100000
	      << "."
	      << BOOST_VERSION / 100 % 1000
	      << "."
	      << BOOST_VERSION % 100
	      << std::endl;

	try {
		boost::archive::binary_iarchive ia(ifs);
		ia >> (*this) ;
	} catch (boost::archive::archive_exception &e) {
		std::cerr << "Argus/ClothFrictionData: cannot read " << file << "; " << e.what() << std::endl ;
		return false ;
	}

	return true ;
}

//! Assemble H matrix et w vector such the relative vel is u = H v + w,
void ClothFrictionData::makeRelativeVelocityAffineMap (bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > & H,
        DynVec &w , bool snapToVertices) const
{

	 const int nContacts = this->nContacts() ;
	 H.reserve( 6 * nContacts ) ;
	 H.setRows( nContacts ) ;
	 H.setCols( M.colsOfBlocks() ) ;

	 w.setZero( 3 * nContacts ) ;

#pragma omp parallel for
	 for( int i = 0 ; i < nContacts ; ++i ) {

		 const Mat3& rotation = contactBasis.diagonal(i).transpose() ;

		 if(snapToVertices) {
			 int k ;
			 coordsA.col(i).maxCoeff(&k);
			 H.insert( i, indicesA(k,i) ) = rotation ;
		 } else {
			 for( int k = 0 ; k< 3 ; ++k ) {
				 H.insert( i, indicesA(k,i) ) = coordsA(k,i) * rotation ;
			 }
		 }

		 w.segment<3>(3*i) = rotation * affineVel.col(i) ;

		 if( indicesB(0,i) != -1 ) {

			 // Assumes no contact between adjacent triangles,
			 //  so we get 6 distinct blocks in the ith row

			 if(snapToVertices) {
				 int k ;
				 coordsB.col(i).maxCoeff(&k);
				 H.insert( i, indicesB(k,i) ) = - rotation ;
			 } else {
				 for( int k = 0 ; k< 3 ; ++k ) {
					 H.insert( i, indicesB(k,i) ) = - coordsB(k,i) * rotation ;
				 }
			 }
		 }
	 }

	 H.finalize() ;
//	 H.prune(1.e-16) ;

}

void ClothFrictionData::applyAdhesion( const DynVec& adh )
{
	assert( adh.rows() == nContacts() ) ;

	//  u = HM-1H' r + HM-1 f
	// (r+adh n,u) in Cmu
	//  rt = r + adh ;  r = rt - adh ; ft = f - H'adh
	//  u = HM-1H' rr + HM-1 ft
	// (rt,u) in Cmu

	// Simply add ( -H' adh e_n )  to f

	bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > H ;
	DynVec w ;
	makeRelativeVelocityAffineMap( H, w );

	w.setZero() ;
	Eigen::Map< DynVec, 0, Eigen::InnerStride<3> >(w.data(), nContacts())
	        = adh ;

	f -= H.transpose() * w ;
}

void ClothFrictionData::alignTangentsWithVelocity(const DynVec &v)
{

	bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > H ;
	DynVec w ;
	makeRelativeVelocityAffineMap( H, w );

	const DynVec u = H*v + w;

	const int n = nContacts();

#pragma omp parallel for
	for( int i = 0 ; i < n ; ++i )
	{
		const Eigen::Vector2d uT = u.segment<2>(3*i+1);

		if( uT.squaredNorm() > 1.e-12 )
		{
			const Mat3& loc2world = contactBasis.diagonal(i).transpose() ;

			const Eigen::Vector3d t1 = (loc2world * Eigen::Vector3d(0, uT[0], uT[1])).normalized();
			contactBasis.diagonal(i).col(1) = t1 ;
			contactBasis.diagonal(i).col(2) = contactBasis.diagonal(i).col(0).cross(t1) ;

		}
	}

}


bool ClothFrictionData::dump(const std::string &file) const
{
	std::ofstream ofs( file.c_str() );

	try {
		boost::archive::binary_oarchive oa(ofs);
		oa << (*this) ;
	} catch (boost::archive::archive_exception &e) {
		std::cerr << "Argus/ClothFrictionData: cannot write " << file << "; " << e.what() << std::endl ;
		return false ;
	}

	return true ;
}




} // argus
