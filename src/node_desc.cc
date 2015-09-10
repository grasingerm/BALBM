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

namespace balbm {

namespace d2q9 {

//! Virtual destructor definition
AbstractNodeDesc::~AbstractNodeDesc() {}

//! D2Q9 base class streaming
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
inline void AbstractNodeDesc::stream(Lattice &lat, const unsigned i,
                                     const unsigned j) const {
#ifdef BALBM_CHECK_BOUNDS_STREAMING
  stream_with_bcheck_(lat, i, j);
#else
  stream_(lat, i, j);
#endif
}

//! D2Q9 base class collide and bound
//!
//! \param
inline void
AbstractNodeDesc::collide_and_bound(Lattice &lat, const CollisionManager &cman,
                                    const unsigned i, const unsigned j) const
    noexcept {
  collide_and_bound_(lat, cman, i, j);
}

//! Virtual destructor definition
AbstractNodeActive::~AbstractNodeActive() {}

//! D2Q9 streaming for a typical active node
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
//! \TODO should we do bounds checking here?
//! \TODO should we unroll this loop?
//! \TODO should this be multithreaded?
void AbstractNodeActive::stream_(Lattice &lat, const unsigned i,
                                 const unsigned j) const {
  const unsigned nk = lat.num_k();
  unsigned i_next, j_next;

  for (unsigned k = 0; k < nk; ++k) {
    j_next = j + lat.c(k, 0);
    i_next = i + lat.c(k, 1);
    assert(lat.in_bounds(i_next, j_next));

    lat.ft_(i_next, j_next, k) = lat.f_(i, j, k);
  }
}

//! D2Q9 streaming for a typical active node with bounds checking
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
//! \TODO should we do bounds checking here?
//! \TODO should we unroll this loop?
//! \TODO should this be multithreaded?
void AbstractNodeActive::stream_with_bcheck_(Lattice &lat, const unsigned i,
                                             const unsigned j) const {
  const unsigned nk = lat.num_k();
  unsigned i_next, j_next;

  for (unsigned k = 0; k < nk; ++k) {
    j_next = j + lat.c(k, 0);
    i_next = i + lat.c(k, 1);

    lat.check_bounds(i_next, j_next);

    lat.ft_(i_next, j_next, k) = lat.f_(i, j, k);
  }
}

//! Default implementation of collide and bound for an active node
//!
//! \param lat Lattice
//! \param cman Collision manager
//! \param mmap Multiscale map
//! \param i Index in the x-direction
//! \param j Index in the y-direction
void AbstractNodeActive::collide_and_bound_(
    Lattice &lat, IncompFlowMultiscaleMap &mmap,
    const IncompFlowCollisionManager &cman, const unsigned i,
    const unsigned j) const {
  cman.collide(lat, mmap, i, j);
}

//! D2Q9 streaming for a west facing node
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
void NodeWestFacingWall::stream_(Lattice &lat, const unsigned i,
                                 const unsigned j) const noexcept {
#ifndef NDEBUG
  const static unsigned stream_directions[] = {2, 3, 4, 6, 7};
  const static unsigned n =
      sizeof(stream_directions) / sizeof(stream_directions[0]);
  unsigned k;

  for (unsigned i = 0; i < n; ++i) {
    k = stream_directions[i];
    j_next = j + lat.c(k, 0);
    i_next = i + lat.c(k, 1);
    assert(lat.in_bounds(i_next, j_next));
  }
#endif

  lat.ft_(i + lat.c(2, 1), j + lat.c(2, 0), 2) = lat.f_(i, j, 2);
  lat.ft_(i + lat.c(3, 1), j + lat.c(3, 0), 3) = lat.f_(i, j, 3);
  lat.ft_(i + lat.c(4, 1), j + lat.c(4, 0), 4) = lat.f_(i, j, 4);
  lat.ft_(i + lat.c(6, 1), j + lat.c(6, 0), 6) = lat.f_(i, j, 6);
  lat.ft_(i + lat.c(7, 1), j + lat.c(7, 0), 7) = lat.f_(i, j, 7);
}

//! D2Q9 streaming for a west facing node with bounds checking
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
void NodeWestFacingWall::stream_with_bcheck_(Lattice &lat, const unsigned i,
                                             const unsigned j) const {
  const static unsigned stream_directions[] = {2, 3, 4, 6, 7};
  const static unsigned n =
      sizeof(stream_directions) / sizeof(stream_directions[0]);
  unsigned k;

  for (unsigned i = 0; i < n; ++i) {
    k = stream_directions[i];
    j_next = j + lat.c(k, 0);
    i_next = i + lat.c(k, 1);

    lat.check_bound(i_next, j_next);
    lat.ft_(i_next, j_next, k) = lat.f_(i, j, k);
  }
}

//! Forward collision call to collision manager then enforce no slip BC
//!
//! \param lat Lattice
//! \param cman Collision manager
//! \param mmap Multiscale map
//! \param i Index in the x-direction
//! \param j Index in the y-direction
void NodeWestFacingWall::collide_and_bound_(
    Lattice &lat, IncompFlowMultiscaleMap &mmap,
    const IncompFlowCollisionManager &cman, const unsigned i,
    const unsigned j) const {
  cman.collide(lat, mmap, i, j);
  lat.f_(i, j, 3) = lat.f_(i, j, 1);
  lat.f_(i, j, 6) = lat.f_(i, j, 8);
  lat.f_(i, j, 7) = lat.f_(i, j, 5);
}

//! D2Q9 streaming for a west facing node
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
void NodeSouthFacingWall::stream_(Lattice &lat, const unsigned i,
                                  const unsigned j) const noexcept {
#ifndef NDEBUG
  const static unsigned stream_directions[] = {1, 3, 4, 7, 8};
  const static unsigned n =
      sizeof(stream_directions) / sizeof(stream_directions[0]);
  unsigned k;

  for (unsigned i = 0; i < n; ++i) {
    k = stream_directions[i];
    j_next = j + lat.c(k, 0);
    i_next = i + lat.c(k, 1);
    assert(lat.in_bounds(i_next, j_next));
  }
#endif

  lat.ft_(i + lat.c(1, 1), j + lat.c(1, 0), 1) = lat.f_(i, j, 1);
  lat.ft_(i + lat.c(3, 1), j + lat.c(3, 0), 3) = lat.f_(i, j, 3);
  lat.ft_(i + lat.c(4, 1), j + lat.c(4, 0), 4) = lat.f_(i, j, 4);
  lat.ft_(i + lat.c(7, 1), j + lat.c(7, 0), 7) = lat.f_(i, j, 7);
  lat.ft_(i + lat.c(8, 1), j + lat.c(8, 0), 8) = lat.f_(i, j, 8);
}

//! D2Q9 streaming for a west facing node with bounds checking
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
void NodeWestFacingWall::stream_with_bcheck_(Lattice &lat, const unsigned i,
                                             const unsigned j) const {
  const static unsigned stream_directions[] = {1, 3, 4, 7, 8};
  const static unsigned n =
      sizeof(stream_directions) / sizeof(stream_directions[0]);
  unsigned k;

  for (unsigned i = 0; i < n; ++i) {
    k = stream_directions[i];
    j_next = j + lat.c(k, 0);
    i_next = i + lat.c(k, 1);

    lat.check_bound(i_next, j_next);
    lat.ft_(i_next, j_next, k) = lat.f_(i, j, k);
  }
}

//! Forward collision call to collision manager then enforce no slip BC
//!
//! \param lat Lattice
//! \param cman Collision manager
//! \param mmap Multiscale map
//! \param i Index in the x-direction
//! \param j Index in the y-direction
void NodeWestFacingWall::collide_and_bound_(
    Lattice &lat, IncompFlowMultiscaleMap &mmap,
    const IncompFlowCollisionManager &cman, const unsigned i,
    const unsigned j) const {
  cman.collide(lat, mmap, i, j);
  lat.f_(i, j, 4) = lat.f_(i, j, 2);
  lat.f_(i, j, 7) = lat.f_(i, j, 5);
  lat.f_(i, j, 8) = lat.f_(i, j, 6);
}

} // namespace d2q9

} // namespace balbm
