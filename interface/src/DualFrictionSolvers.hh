#ifndef ARGUS_DUAL_FRICTION_SOLVERS_HH
#define ARGUS_DUAL_FRICTION_SOLVERS_HH

#include "MatrixTypes.hh"

#include <bogus/Core/Block.hpp>

namespace argus {

struct ClothFrictionData ;
struct SolverOptions ;
struct SolverStats ;

class DualFrictionSolver {
public:

	typedef bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > HType ;

	DualFrictionSolver(
	        const ClothFrictionData& data,
	        const HType &H,
	        const DynVec &w
	        ) ;

	~DualFrictionSolver();

	double solve( const SolverOptions& options, DynVec &v, DynVec &r, SolverStats& stats ) const ;

private:

	void factorize() ;

	const ClothFrictionData& m_data ;
	const HType&  m_H ;
	const DynVec& m_w ;

	struct MFactorization ;
	const MFactorization * m_fac ;
};

class ICAFrictionSolver {
public:

	typedef bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > HType ;
	typedef bogus::SparseBlockMatrix< Mat3, bogus::SYMMETRIC > DType ;

	ICAFrictionSolver(
	        const ClothFrictionData& data,
	        const HType &H,
	        const DynVec &w
	        ) ;

	double solve( const SolverOptions& options, DynVec &v, DynVec &r, SolverStats& stats ) const ;

private:


	const ClothFrictionData& m_data ;
	const HType&  m_H ;
	const DynVec& m_w ;
};


} //argus


#endif
