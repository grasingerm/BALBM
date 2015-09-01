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
// this program.  If not, see <http://www.gnu.org/licenses/>

#include "multiscale_map.hh"

namespace balbm
{

namespace d2q9
{

//! Virtual destructor for base class
AbstractMultiscaleMap::~AbstractMultiscaleMap() {}

//! Map particle distribution functions to density
//!
//! \param lat D2Q9 lattice
AbstractMultiscaleMap::map_to_macro_(const Lattice& lat)
{
  for (unsigned j = 0; j < num_y(); ++j)
    for (unsigned i = 0; i < num_x(); ++i)
      map_to_macro_(i, j);
}

//! Map particle distribution functions to density
//!
//! \param lat D2Q9 lattice
//! \param i x-coord of node
//! \param j y-coord of node
void AbstractMultiscaleMap::map_to_macro_
  (const Lattice& lat, const unsigned i, const unsigned j)
{
  const nk = lat.num_k();
  rho_(i, j) = 0.;
  for (unsigned k = 0; i < nk; ++k)
    rho_(i, j) += lat.f(k, i, j);
}

//! Map particle distribution functions to incompressible flow macroscopic 
//! variables
//!
//! \param lat D2Q9 lattice
//! \param i x-coord of node
//! \param j y-coord of node
void IncompFlowMultiscaleMap::_map_to_macro
  (const Lattice& lat, const unsigned i, const unsigned j)
{
  const nk = lat.num_k();
  rho_(i, j) = 0.;
  u_(i, j, 0) = 0.;
  u_(i, j, 1) = 0.;

  for (unsigned k = 0; k < nk; ++k)
  {
    rho_(i, j) += lat.f(i, j, k);
    u_(0, i, j) += lat.k(0, k) * lat.f(k, i, j);
    u_(1, i, j) += lat.k(1, k) * lat.f(k, i, j);
  }

  u_(0, i, j) /= rho_(i, j);
  u_(1, i, j) /= rho_(i, j);
}

} // namespace d2q9

} // namespace balbm
