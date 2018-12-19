#ifndef ARGUS_NODAL_CONTACT_HH
#define ARGUS_NODAL_CONTACT_HH

#include "MatrixTypes.hh"

namespace argus {

struct ClothFrictionData ;
struct SolverOptions ;
struct SolverStats ;

class NodalContact {
public:

	NodalContact( const ClothFrictionData& data )
		: m_data(data)
	{}


	double solve(const SolverOptions& options, DynVec &v, DynVec &r, SolverStats& stats ) const ;

	double solve(const SolverOptions& options, DynVec &v, SolverStats& stats ) const ;

private:

	const ClothFrictionData& m_data ;
};


} //argus

#endif
