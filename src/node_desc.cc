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

#include "collision_manager.hh"
#include "lattice.hh"
#include "node_desc.hh"
#include <cassert>

namespace balbm {

namespace d2q9 {

//! Virtual destructor definition
AbstractNodeDesc::~AbstractNodeDesc() {}

//! D2Q9 base class streaming
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
void AbstractNodeDesc::stream(Lattice &lat, const unsigned i,
                              const unsigned j) const {
#ifdef BALBM_CHECK_BOUNDS_STREAMING
  stream_with_bcheck_(lat, i, j);
#else
  stream_(lat, i, j);
#endif
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
                                 const unsigned j) const noexcept {
  const unsigned nk = lat.num_k();
  unsigned i_next, j_next;

  for (unsigned k = 0; k < nk; ++k) {
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);
    assert(lat.in_bounds(i_next, j_next));

    lat.ft(i_next, j_next, k) = lat.f(i, j, k);
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
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);

    lat.check_bounds(i_next, j_next);

    lat.ft(i_next, j_next, k) = lat.f(i, j, k);
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
  unsigned k, i_next, j_next;

  for (unsigned idx = 0; idx < n; ++idx) {
    k = stream_directions[idx];
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);
    assert(lat.in_bounds(i_next, j_next));
  }
#endif

  lat.ft(i + lat.c(2, 0), j + lat.c(2, 1), 2) = lat.f(i, j, 2);
  lat.ft(i + lat.c(3, 0), j + lat.c(3, 1), 3) = lat.f(i, j, 3);
  lat.ft(i + lat.c(4, 0), j + lat.c(4, 1), 4) = lat.f(i, j, 4);
  lat.ft(i + lat.c(6, 0), j + lat.c(6, 1), 6) = lat.f(i, j, 6);
  lat.ft(i + lat.c(7, 0), j + lat.c(7, 1), 7) = lat.f(i, j, 7);
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
  unsigned k, i_next, j_next;

  for (unsigned idx = 0; idx < n; ++idx) {
    k = stream_directions[idx];
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);

    lat.check_bounds(i_next, j_next);
    lat.ft(i_next, j_next, k) = lat.f(i, j, k);
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
  lat.f(i, j, 3) = lat.f(i, j, 1);
  lat.f(i, j, 6) = lat.f(i, j, 8);
  lat.f(i, j, 7) = lat.f(i, j, 5);
}

//! D2Q9 streaming for a south facing node
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
  unsigned k, i_next, j_next;

  for (unsigned idx = 0; idx < n; ++idx) {
    k = stream_directions[idx];
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);
    assert(lat.in_bounds(i_next, j_next));
  }
#endif

  lat.ft(i + lat.c(1, 0), j + lat.c(1, 1), 1) = lat.f(i, j, 1);
  lat.ft(i + lat.c(3, 0), j + lat.c(3, 1), 3) = lat.f(i, j, 3);
  lat.ft(i + lat.c(4, 0), j + lat.c(4, 1), 4) = lat.f(i, j, 4);
  lat.ft(i + lat.c(7, 0), j + lat.c(7, 1), 7) = lat.f(i, j, 7);
  lat.ft(i + lat.c(8, 0), j + lat.c(8, 1), 8) = lat.f(i, j, 8);
}

//! D2Q9 streaming for a south facing node with bounds checking
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
void NodeSouthFacingWall::stream_with_bcheck_(Lattice &lat, const unsigned i,
                                              const unsigned j) const {
  const static unsigned stream_directions[] = {1, 3, 4, 7, 8};
  const static unsigned n =
      sizeof(stream_directions) / sizeof(stream_directions[0]);
  unsigned k, i_next, j_next;

  for (unsigned idx = 0; idx < n; ++idx) {
    k = stream_directions[idx];
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);

    lat.check_bounds(i_next, j_next);
    lat.ft(i_next, j_next, k) = lat.f(i, j, k);
  }
}

//! Forward collision call to collision manager then enforce no slip BC
//!
//! \param lat Lattice
//! \param cman Collision manager
//! \param mmap Multiscale map
//! \param i Index in the x-direction
//! \param j Index in the y-direction
void NodeSouthFacingWall::collide_and_bound_(
    Lattice &lat, IncompFlowMultiscaleMap &mmap,
    const IncompFlowCollisionManager &cman, const unsigned i,
    const unsigned j) const {
  cman.collide(lat, mmap, i, j);
  lat.f(i, j, 4) = lat.f(i, j, 2);
  lat.f(i, j, 7) = lat.f(i, j, 5);
  lat.f(i, j, 8) = lat.f(i, j, 6);
}

//! D2Q9 streaming for a east facing node
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
void NodeEastFacingWall::stream_(Lattice &lat, const unsigned i,
                                 const unsigned j) const noexcept {
#ifndef NDEBUG
  const static unsigned stream_directions[] = {1, 2, 4, 5, 8};
  const static unsigned n =
      sizeof(stream_directions) / sizeof(stream_directions[0]);
  unsigned k, i_next, j_next;

  for (unsigned idx = 0; idx < n; ++idx) {
    k = stream_directions[idx];
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);
    assert(lat.in_bounds(i_next, j_next));
  }
#endif

  lat.ft(i + lat.c(1, 0), j + lat.c(1, 1), 1) = lat.f(i, j, 1);
  lat.ft(i + lat.c(2, 0), j + lat.c(2, 1), 2) = lat.f(i, j, 2);
  lat.ft(i + lat.c(4, 0), j + lat.c(4, 1), 4) = lat.f(i, j, 4);
  lat.ft(i + lat.c(5, 0), j + lat.c(5, 1), 5) = lat.f(i, j, 5);
  lat.ft(i + lat.c(8, 0), j + lat.c(8, 1), 8) = lat.f(i, j, 8);
}

//! D2Q9 streaming for a east facing node with bounds checking
//!
//! \param lat D2Q9 lattice
//! \param i y-coord of node
//! \param j x-coord of node
void NodeEastFacingWall::stream_with_bcheck_(Lattice &lat, const unsigned i,
                                             const unsigned j) const {
  const static unsigned stream_directions[] = {1, 2, 4, 5, 8};
  const static unsigned n =
      sizeof(stream_directions) / sizeof(stream_directions[0]);
  unsigned k, i_next, j_next;

  for (unsigned idx = 0; idx < n; ++idx) {
    k = stream_directions[idx];
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);

    lat.check_bounds(i_next, j_next);
    lat.ft(i_next, j_next, k) = lat.f(i, j, k);
  }
}

//! Forward collision call to collision manager then enforce no slip BC
//!
//! \param lat Lattice
//! \param cman Collision manager
//! \param mmap Multiscale map
//! \param i Index in the x-direction
//! \param j Index in the y-direction
void NodeEastFacingWall::collide_and_bound_(
    Lattice &lat, IncompFlowMultiscaleMap &mmap,
    const IncompFlowCollisionManager &cman, const unsigned i,
    const unsigned j) const {
  cman.collide(lat, mmap, i, j);
  lat.f(i, j, 1) = lat.f(i, j, 3);
  lat.f(i, j, 5) = lat.f(i, j, 7);
  lat.f(i, j, 8) = lat.f(i, j, 6);
}

// TODO: prob micro-opt, BUT what-if we skip all 0 for streaming???

//! D2Q9 streaming for a north facing node
//!
//! \param lat D2Q9 lattice
//! \param i x-coord of node
//! \param j y-coord of node
void NodeNorthFacingWall::stream_(Lattice &lat, const unsigned i,
                                  const unsigned j) const noexcept {
#ifndef NDEBUG
  const static unsigned stream_directions[] = {1, 3, 4, 7, 8};
  const static unsigned n =
      sizeof(stream_directions) / sizeof(stream_directions[0]);
  unsigned k, i_next, j_next;

  for (unsigned idx = 0; idx < n; ++idx) {
    k = stream_directions[idx];
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);
    assert(lat.in_bounds(i_next, j_next));
  }
#endif

  lat.ft(i + lat.c(1, 0), j + lat.c(1, 1), 1) = lat.f(i, j, 1);
  lat.ft(i + lat.c(2, 0), j + lat.c(2, 1), 2) = lat.f(i, j, 2);
  lat.ft(i + lat.c(3, 0), j + lat.c(3, 1), 3) = lat.f(i, j, 3);
  lat.ft(i + lat.c(5, 0), j + lat.c(5, 1), 5) = lat.f(i, j, 5);
  lat.ft(i + lat.c(6, 0), j + lat.c(6, 1), 6) = lat.f(i, j, 6);
}

//! D2Q9 streaming for a north facing node with bounds checking
//!
//! \param lat D2Q9 lattice
//! \param i x-coord of node
//! \param j y-coord of node
void NodeNorthFacingWall::stream_with_bcheck_(Lattice &lat, const unsigned i,
                                              const unsigned j) const {
  const static unsigned stream_directions[] = {1, 2, 3, 5, 6};
  const static unsigned n =
      sizeof(stream_directions) / sizeof(stream_directions[0]);
  unsigned k, i_next, j_next;

  for (unsigned idx = 0; idx < n; ++idx) {
    k = stream_directions[idx];
    i_next = i + lat.c(k, 0);
    j_next = j + lat.c(k, 1);

    lat.check_bounds(i_next, j_next);
    lat.ft(i_next, j_next, k) = lat.f(i, j, k);
  }
}

//! Forward collision call to collision manager then enforce no slip BC
//!
//! \param lat Lattice
//! \param cman Collision manager
//! \param mmap Multiscale map
//! \param i Index in the x-direction
//! \param j Index in the y-direction
void NodeNorthFacingWall::collide_and_bound_(
    Lattice &lat, IncompFlowMultiscaleMap &mmap,
    const IncompFlowCollisionManager &cman, const unsigned i,
    const unsigned j) const {
  cman.collide(lat, mmap, i, j);
  lat.f(i, j, 2) = lat.f(i, j, 4);
  lat.f(i, j, 5) = lat.f(i, j, 7);
  lat.f(i, j, 6) = lat.f(i, j, 8);
}

//! Constructor for periodic boundary condition node
//!
//! \param i_next Index in x-direction of node to copy particle distributions to
//! \param j_next Index in y-direction of node to copy particle distributions to
//! \param ks Array of indexes for lattice directions to be copied
//! \param nk Number of elements in ks
//! \return Periodic boundary condition node descriptor
NodePeriodic::NodePeriodic(const unsigned i_next, const unsigned j_next,
                           const unsigned *ks, const unsigned nk)
    : i_next_(i_next), j_next_(j_next),
      ks_(std::unique_ptr<unsigned[]>(new unsigned[nk])), nk_(nk) {
  assert(static_cast<int>(nk) >= 0 &&
         nk <= 9); // TODO: consider turning this into a throw
  for (unsigned k = 0; k < nk; ++k) {
    assert(ks[k] < 9 && static_cast<int>(ks[k]) >= 0);
    ks_[k] = ks[k];
  }
}

//! Forward collision call to collision manager then enforce periodic BC
//!
//! \param lat Lattice
//! \param cman Collision manager
//! \param mmap Multiscale map
//! \param i Index in the x-direction
//! \param j Index in the y-direction
void NodePeriodic::collide_and_bound_(Lattice &lat,
                                      IncompFlowMultiscaleMap &mmap,
                                      const IncompFlowCollisionManager &cman,
                                      const unsigned i,
                                      const unsigned j) const {
  cman.collide(lat, mmap, i, j);
  assert(lat.in_bounds(i_next_, j_next_)); // TODO: throw?
  for (unsigned k = 0; k < nk_; ++k)
    lat.f(i_next_, j_next_, k) = lat.f(i, j, k);
}

} // namespace d2q9

} // namespace balbm
