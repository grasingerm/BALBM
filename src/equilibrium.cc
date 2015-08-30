#include "equilibrium.hh"

namespace balbm
{

namespace d2q9
{

//! Virtual destructor definition
AbstractEqFunct::~AbstractEqFunct() {}

//! Equilibrium distribution function for incompressible flow
//!
//! \param
double IncompFlowEqFunct::_f const(const Lattice& lat, 
                                   const AbstractMultiscaleMap& abmmap,
                                   const unsigned k)
{
  auto inmmap = dynamic_cast<const IncompFlowMultiscaleMap&>(abmmap);
  double k_dot_u = lat.k(k, 0) * inmmap.u
}

}

}
