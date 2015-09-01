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

#include "node_desc.hh"
#include <cassert>
#include <exception>
#include <sstream>

namespace balbm
{

namespace d2q9
{

//! Virtual destructor definition
AbstractNodeDesc::~AbstractNodeDesc() {}

//! D2Q9 base class streaming
//!
//! \param lat D2Q9 lattice
//! \param i x-coord of node
//! \param j y-coord of node
inline void AbstractNodeDesc::stream (Lattice& lat, const unsigned i, const unsigned j) 
  const
{
#ifdef BALBM_CHECK_BOUNDS_STREAMING
  _stream_with_bcheck(lat, i, j);
#else
  _stream(lat, i, j);
#endif
}

//! D2Q9 base class collide and bound
//!
//! \param
inline void AbstractNodeDesc::collide_and_bound
  (Lattice& lat, const CollisionManager& cman, const unsigned i, 
   const unsigned j) const noexcept
{
  collide_and_bound_(lat, cman, i, j);
}

//! Virtual destructor definition
AbstractNodeActive::~AbstractNodeActive() {}

//! D2Q9 streaming for a typical active node
//!
//! \param lat D2Q9 lattice 
//! \param i x-coord of node
//! \param j y-coord of node
//! \TODO should we do bounds checking here?
//! \TODO should we unroll this loop?
//! \TODO should this be multithreaded?
void AbstractNodeActive::stream_(Lattice& lat, const unsigned i, 
                                 const unsigned j) const
{
  const unsigned nk = lat.num_k();
  unsigned i_next, j_next;

  for (unsigned k = 0; k < nk; ++k)
  {
    i_next = i + lat.k(k, 0);
    j_next = j + lat.k(k, 1);
    assert(i_next < lat.num_x() && i_next >= 0);
    assert(j_next < lat.num_y() && j_next >= 0);

    lat.ft_(i_next, j_next, k) = lat.f_(i, j, k);
  }
}

//! D2Q9 streaming for a typical active node with bounds checking
//!
//! \param lat D2Q9 lattice 
//! \param i x-coord of node
//! \param j y-coord of node
//! \TODO should we do bounds checking here?
//! \TODO should we unroll this loop?
//! \TODO should this be multithreaded?
void AbstractNodeActive::stream_with_bcheck_(Lattice& lat, const unsigned i, 
                                             const unsigned j) const
{
  const unsigned nk = lat.num_k();
  unsigned i_next, j_next;

  for (unsigned k = 0; k < nk; ++k)
  {
    i_next = i + lat.k(k, 0);
    j_next = j + lat.k(k, 1);

    // check that streaming stays in bounds
    if (!(i_next < lat.num_x() && i_next >= 0 
          && j_next < lat.num_y() && j_next >= 0))
    {
      std::ostringstream oss;
      oss << "Ill-defined boundaries. Particles streamed out of bounds from "
          << "node (" << i << " ," << j << ") to (" << i_next << ", "
          << j_next << "). Check boundary conditions.";
      throw std::out_of_range (oss.str());
    }

    lat.ft_(i_next, j_next, k) = lat.f_(i, j, k);
  }
}

void AbstractNodeActive::_collide_and_bound
  (Lattice& lat, const CollisionManager& cman, const MultiscaleMap& mmap, 
   const unsigned i, const unsigned j) const
{
  // code for collisions
}

}
