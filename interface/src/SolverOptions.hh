#ifndef ARGUS_SOLVER_OPTIONS_HH
#define ARGUS_SOLVER_OPTIONS_HH

#include <limits>
#include <vector>

namespace argus {

//! Tolerance and solver tuning parameters
struct SolverOptions {
	enum Algorithm {
		NodalContact = 0,              //!< SCA poster algorithm on primal formulation
		AlternatingMinimization = 1,   //!< bogus DualAMA primal-dual algorithm
		DualProjectedGradient = 2,     //!< PG on dual formulation -- req. matrix factorization
		DualGaussSeidel = 3,            //!< GS on dual formulation -- req. matrix factorization
		NodalSelfContact = 4,     //!< testing out nodal self intersections //
		TotalContact = 5,			// Test handling entire impacts
		ICAGaussSeidel = 6             //!< GS with Iterative constraint anticipation
	};

	Algorithm algorithm ;
	double    tolerance ;
	bool useInfinityNorm ;

	unsigned maxIterations ;	  //Inner
	unsigned maxOuterIterations ; //For Cadoux algorithm
	unsigned lineSearchIterations ; // For AMA and dome PG variants

	int projectedGradientVariant ;

	double amaFpStepSize ;          // AMA fixed-point step size
	double amaProjStepSize ;        // AMA projection step size


	double nodalConstraintStepSize ;
	int    nodalConstraintIterations ;
	double nodalGSAutoRegularization ;
	double nodalGSLocalTolerance ;

	bool   faceted ;

	double AScale;
	double bScale;
	double evalScale;

	std::vector<double> m_inv;

	bool   useColoring;

	SolverOptions()
	    : algorithm(NodalContact),
	      tolerance( 1.e-8 ),
	      useInfinityNorm( false ),
	      maxIterations( 500 ),
	      maxOuterIterations( 25 ),
	      lineSearchIterations( 8 ),
	      projectedGradientVariant( -1 ),
	      amaFpStepSize(.05),
	      amaProjStepSize(1.e-3),
	      nodalConstraintStepSize(1.5),
	      nodalConstraintIterations(0),
	      nodalGSAutoRegularization(1.e-3),
	      nodalGSLocalTolerance(
	          std::pow( std::numeric_limits< double >::epsilon(), .75 ) ),
	      faceted(false),
	      AScale(1),
	      bScale(1),
	      evalScale(1),
	      useColoring(false)

	{}
};

} //argus

#endif
