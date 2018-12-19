#ifndef ARGUS_CLOTH_FRICTION_SOLVER_HH
#define ARGUS_CLOTH_FRICTION_SOLVER_HH

#include "MatrixTypes.hh"

#include "SolverOptions.hh"
#include "SolverStats.hh"

namespace argus {

struct ClothFrictionData ;

//! Class for solving a ClothFrictionData contact problem
class ClothFrictionSolver {
public:

	ClothFrictionSolver( const ClothFrictionData& data )
	    : m_data(data)
	{}

	//! Solves the friction problem m_data
	/*!
	 * \param options [in] Solver configuration: tolerance, choice of algorithm, etc
	 * \param v       [in/out] Velocity variable: vector of size m_data.M.rows(), properly initialized (eg 0 or last step velocities)
	 * \param stats   [out] Info about the solver convergence
	 * \return the solver minimization error, or -1 if an error occured
	 */
	double solve(const SolverOptions& options, DynVec& v, SolverStats &stats) const ;
	double solve(const SolverOptions& options, DynVec& v, DynVec& r, SolverStats &stats) const ;

public:
	void linearSolve(DynVec& v) const ;

	double contactSolve(const SolverOptions& options, DynVec& v, DynVec& r, SolverStats &stats) const ;

	double solveContact(const SolverOptions& options, DynVec& v, DynVec& r, SolverStats &stats) const ;

	const ClothFrictionData& m_data ;

};


} //argus

#endif
