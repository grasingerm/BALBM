// Complex flow simulator using lattice Boltzmann method
// Copyright (C) 2015 Matthew Grasinger
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// A copy of the GNU General Public License is at the root directory of

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
double IncompFlowEqFunct::f_ const(const Lattice& lat, 
                                   const AbstractMultiscaleMap& abmmap,
                                   const unsigned k)
{
  auto inmmap = dynamic_cast<const IncompFlowMultiscaleMap&>(abmmap);
  double k_dot_u = lat.k(k, 0) * inmmap.u
}

}

}
