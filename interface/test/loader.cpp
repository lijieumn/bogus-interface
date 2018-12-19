#include "ClothFrictionData.hh"
#include "ClothFrictionSolver.hh"

#include <bogus/Core/Block.io.hpp>
#include <bogus/Core/Block.impl.hpp>

#include <fstream>

static void usage( const char *name )
{
	std::cout << "Usage: " << name
	          << " problem_file [options] "
	          << "\nOffline cloth friction problem solver. "
	          << "\n\n" ;

	std::cout << "Options:\n"
	          << "-h \t Display this help message and exit\n"
	          << "-t tol \t Tolerance \n"
	          << "-i useInfNorm \t Whether to use infinity norm instead of l_2\n"
	          << "-a algo \t Algorithm id in [0,3]  \n"
	          << "-g variant \t PG variant id in [0,4]  \n"
	          << "-n nIters \t Max number of inner solver iterations  \n"
	          << "-o nIters \t Max number of outer fixed-point iterations  \n"
	          << "-l nIters \t Max number of line-search iterations  \n"
	          << "-f size \t AMA fixed-point step size\n"
	          << "-p size \t AMA projection step size\n"
	          << "-c size \t Nodal constraint step size\n"
	          << "-d bool \t Enforce determinicity\n"
	          << "-F bool \t Faceted (for ICA) \n"
	          << "-V path \t Ascii initial velocity file \n"
	          << std::endl ;
}


static void solveSample( const argus::SolverOptions& options,
                         argus::ClothFrictionData& data )
{

	argus::SolverStats   stats   ;

	argus::DynVec v( data.M.rows() ) ;
	v.setZero() ;

	if(options.algorithm == argus::SolverOptions::NodalContact ||
	        options.algorithm == argus::SolverOptions::NodalSelfContact )
		data.findAndDuplicate() ;
	///		std::cout << data.M << std::endl ;
	//		std::cout << data.duplicatedNodes[0] << std::endl ;
	//		std::cout << data.duplicatedNodes[1] << std::endl ;
	//		std::cout << data.f.transpose() << std::endl ;

	argus::ClothFrictionSolver solver(data) ;
	solver.solve( options, v, stats ) ;

	std::cout << "v is: \n" << v.transpose() << std::endl;

	//		argus::DynVec merged ;
	//		data.gatherVelocities( v, merged );
	//		std::cout << "merged as: \n" << merged.transpose() << std::endl;
}


namespace argus
{

struct ClothFrictionEvalFunctor : public EvalFunctor
{

	ClothFrictionEvalFunctor(
	        const SolverOptions &options,
	        const ClothFrictionData& data, DynVec& v, DynVec& r)
	    : m_options(options)
	    , m_data(data)
	    , m_M(data.M)
	    , m_f(data.f)
	    , m_mu(data.mu)
	    , m_v(v)
	    , m_r(r)
	{
		data.makeRelativeVelocityAffineMap( m_H, m_w, true );
	}

	double operator() () const override
	{
		// Typical r vs typical u (~mass)
		constexpr double rho = 1.e-2 ;

		DynVec vv ;

		if( m_v.rows() == m_f.rows())
		{
			vv = m_v ;
		}  else {
			m_data.gatherVelocities( m_v, vv ) ;
		}

		DynVec u = m_H*vv + m_w ;

		double friction_error = 0 ;
		{
			const int n = m_mu.rows();
			double err = 0 ;
#pragma omp parallel for reduction( + : err )
			for( int i = 0 ; i < n ; ++i ) {
				if( m_mu[i] >= 0 )
				{
					// Pi_R+ ( r_n - rho u_N )
					// Pi_B(mu N) ( r_T - rho u_T )

					const Eigen::Vector3d lx = m_r.segment<3>(3*i) ;
					Eigen::Vector3d  ac = lx - rho * u.segment<3>(3*i) ;
					// Normal part
					ac[0] = std::max(0., ac[0])  ;
					//Tangential part
					const double nT = ac.segment<2>(1).norm() ;
					const double rN = std::max(0., lx[0])*m_mu[i];
					if( nT > rN ) {
						ac.segment<2>(1) *= rN/nT ;
					}
					//Error
					ac -= lx ;
					const double lerr = ac.squaredNorm() ;
					err += lerr ;
				}
			}
			friction_error = err / (1+n) ;
		}


		DynVec force_resid = m_M*vv - m_f - m_H.transpose()*m_r ;
		double force_error = force_resid.squaredNorm() / (1+force_resid.rows());

		return  force_error + friction_error ;
	}

	const SolverOptions& m_options;
	const ClothFrictionData& m_data ;

	ClothFrictionData::StiffnessMatrixType m_M ;
	bogus::SparseBlockMatrix< Mat3, bogus::UNCOMPRESSED > m_H ;
	DynVec m_w ;
	DynVec m_f ;
	DynVec m_mu ;

	DynVec& m_v;
	DynVec& m_r;
};
}

namespace
{
void readVelocity( argus::DynVec&v, const char* fname )
{
	std::ifstream in(fname) ;

	int i = 0;
	while( i < v.rows() && in >> v[i++] ) ;

	std::cout << "Read " << i << " velocity components" << std::endl ;
}
}


int main( int argc, const char* argv[] )
{
	argus::SolverOptions options ;
	argus::SolverStats   stats   ;

	const char* problem = 0 ;
	const char* velocityFile = 0 ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
			case 'h':
				usage(argv[0]) ;
				return 0 ;
			case 'a':
				if( ++i == argc ) break ;
				options.algorithm = (argus::SolverOptions::Algorithm) std::atoi( argv[i] ) ;
				break ;
			case 'i':
				if( ++i == argc ) break ;
				options.useInfinityNorm = (bool) std::atoi( argv[i] ) ;
				break ;
			case 'n':
				if( ++i == argc ) break ;
				options.maxIterations = std::atoi( argv[i] ) ;
				break ;
			case 'o':
				if( ++i == argc ) break ;
				options.maxOuterIterations = std::atoi( argv[i] ) ;
				break ;
			case 'l':
				if( ++i == argc ) break ;
				options.lineSearchIterations = std::atoi( argv[i] ) ;
				break ;
			case 'g':
				if( ++i == argc ) break ;
				options.projectedGradientVariant = std::atoi( argv[i] ) ;
				break ;
			case 't':
				if( ++i == argc ) break ;
				options.tolerance = std::strtod( argv[i], 0 ) ;
				break ;
			case 'T':
				if( ++i == argc ) break ;
				stats.exitTolerance = std::strtod( argv[i], 0 ) ;
				break ;
			case 'M':
				if( ++i == argc ) break ;
				stats.exitTime = std::strtod( argv[i], 0 ) ;
				break ;
			case 'f':
				if( ++i == argc ) break ;
				options.amaFpStepSize = std::strtod( argv[i], 0 ) ;
				break ;
			case 'p':
				if( ++i == argc ) break ;
				options.amaProjStepSize = std::strtod( argv[i], 0 ) ;
				break ;
			case 'c':
				if( ++i == argc ) break ;
				options.nodalConstraintStepSize = std::strtod( argv[i], 0 ) ;
				break ;
			case 'd':
				if( ++i == argc ) break ;
				options.useColoring = std::atoi( argv[i] ) ;
				break ;
			case 'F':
				if( ++i == argc ) break ;
				options.faceted = std::atoi( argv[i] ) ;
				break ;
			case 'V':
				if( ++i == argc ) break ;
				velocityFile = argv[i] ;
				break ;

			}
		} else {
			problem = argv[i] ;
		}
	}

	if( !problem ) {

		// If no problem file was specified, then we set up a basic scenario where a vertex from a triangle contacts an external object

		bool ICA = true ;

		argus::SolverOptions options ;
		options.maxOuterIterations = ICA ? 10 : 0 ;
		options.tolerance = 1.e-12 ;
		options.algorithm = ICA
		        ? argus::SolverOptions::ICAGaussSeidel
		        : argus::SolverOptions::NodalSelfContact;
		options.nodalConstraintStepSize   = 1.25;

		options.faceted = true && ICA ;

		{
			argus::ClothFrictionData data;
			data.loadSimpleTest();
			data.applyAdhesion( 10*argus::DynVec::Ones(data.nContacts()) );
			solveSample( options, data ) ;
		}
		{
			argus::ClothFrictionData data;
			data.loadPinTest();
			solveSample( options, data ) ;
		}
		{
			argus::ClothFrictionData data;
			data.loadTwoContactsTest();
			solveSample( options, data ) ;
		}
		{
			argus::ClothFrictionData data;
			data.loadTwoLayersTest();
			solveSample( options, data ) ;
		}
		{
			argus::ClothFrictionData data;
			data.loadSlidingTest();
			solveSample( options, data ) ;
		}

		usage(argv[0]) ;
		return 1 ;
	}

	argus::ClothFrictionData data ;
	if( data.load(problem) ) {

		options.faceted = options.faceted && options.algorithm == argus::SolverOptions::ICAGaussSeidel;

		std::cout << "Number of vertices: " << data.nVertices() << std::endl ;
		std::cout << "Number of contacts: " << data.nContacts() << std::endl ;

		argus::DynVec v, r ;

		stats.functor.reset( new argus::ClothFrictionEvalFunctor(options, data, v, r) );

		if( data.duplicatedNodes.empty() && (
		            options.algorithm == argus::SolverOptions::NodalContact ||
		            options.algorithm == argus::SolverOptions::NodalSelfContact ||
		            options.algorithm == argus::SolverOptions::TotalContact ) )
		{
			data.findAndDuplicate();
		}
		std::cout << "Number of duplicated nodes: " << data.nDuplicatedVertices() << std::endl ;

//		for( unsigned i = 0 ; i < data.duplicatedNodes.size() ; ++i ) {
//			std::cerr <<"(";
//			for( unsigned j = 0 ; j < data.duplicatedNodes[i].size() ; ++j )
//				std::cerr << data.duplicatedNodes[i][j] << " " ;
//			std::cerr <<") ";
//		}
//		std::cerr << std::endl ;


		argus::ClothFrictionSolver solver(data) ;

		v.setZero(data.M.rows()) ;
		r.setZero(data.nContacts() * 3) ;

		if( velocityFile )
			readVelocity(v, velocityFile);

		if(options.faceted)
			data.alignTangentsWithVelocity(v);

		try
		{
			solver.solve( options, v, r, stats ) ;
		} catch (argus::SolverStats::TimeExceeded& ) {
			std::cerr << "Time exceeded, " << "err = " << stats.error << " \t time = " << stats.time << std::endl;
			return 1;
		} catch (argus::SolverStats::ToleranceReached& ) {
			std::cout << "Tolerance reached, " << "err = " << stats.error << " \t time = " << stats.time << std::endl;
		}

//		std::cerr << "Objective error " << (*stats.functor)() << std::endl ;

//		std::cerr << "Velocities: \n" << v.transpose() << std::endl ;

		return 0 ;
	}

	return 1 ;

}
