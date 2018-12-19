#ifndef ARGUS_SOLVER_STATS_HH
#define ARGUS_SOLVER_STATS_HH

#include <bogus/Core/Utils/Timer.hpp>

#include <iostream>
#include <memory>

namespace argus {

struct EvalFunctor
{
	~EvalFunctor() {}
	virtual double operator()() const = 0;
};


//! Info about the error, number of iterations and time taken bu the solver
struct SolverStats
{
	double time ;
	double error ;
	int nIterations ;

	std::unique_ptr<EvalFunctor> functor ;

	struct TimeExceeded {};
	struct ToleranceReached {};
	double exitTolerance ;
	double exitTime ;


	//! Callback from solver algorithms
	void ack( unsigned iter, double err )
	{
#ifndef SILENCE_ARGUS
		std::cout << "ClothFriction iter=" << iter << "\t err =" << err ;
#endif
		error = err ;

		if( functor )
		{
			m_errorEvalTimer.reset();

			error = (*functor)();
#ifndef SILENCE_ARGUS
			std::cout << "\t func = " <<  error ;
#endif

			m_errorEvalTime += m_errorEvalTimer.elapsed();
		}

		nIterations = iter ;
		time = m_timer.elapsed() - m_errorEvalTime ;

#ifndef SILENCE_ARGUS
		std::cout << " \t T= " << time << "\n" ;
#endif

		if( exitTime > 0 && time > exitTime )
			throw TimeExceeded() ;
		if( exitTolerance > 0 && error < exitTolerance )
			throw ToleranceReached() ;
	}

	//! Resets time and counters
	void reset() {
		m_timer.reset();
		time = 0 ;
		nIterations = 0 ;
		error = -1 ;
		m_errorEvalTime = 0;
	}

	SolverStats()
	    : exitTolerance(0), exitTime(0)
	{
		reset() ;
	}

private:
	bogus::Timer m_timer ;
	bogus::Timer m_errorEvalTimer ;
	double m_errorEvalTime;
};

} //argus

#endif
