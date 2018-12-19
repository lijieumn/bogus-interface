#include "RevCoulombLaw.impl.hh"

using namespace bogus ;

namespace argus {

template class RevCoulombLaw< 3u, double, local_soc_solver::RevHybrid > ;
template class RevCoulombLaw< 3u, double, local_soc_solver::PureNewton > ;
template class RevCoulombLaw< 3u, double, local_soc_solver::Hybrid > ;
template class RevCoulombLaw< 3u, double, local_soc_solver::PureEnumerative > ;

} //argus

