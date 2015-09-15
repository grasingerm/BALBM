#ifndef NODE_DESC_HH
#define NODE_DESC_HH

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

// Q: why wrap collide and bound into a single virtual function?
// A: because whether or not a collision occurs at a node is part of its BC.
//    AND because this only requires one vtable lookup instead of two

#include "balbm_config.hh"
#include <algorithm>
#include <memory>

namespace balbm {

namespace d2q9 {

class Lattice;
class IncompFlowMultiscaleMap;
class IncompFlowCollisionManager;

//! \class AbstractNodeDesc
//!
//! \brief Abstract base class for  node descriptors
//!
//! Provide polymorphic behavior for each node in the lattice based on its
//! physical "status" as a node, e.g. change behavior of streaming and collision
//! steps in order to simulate appropriate physics and boundary conditions
class AbstractNodeDesc {
public:
  void stream(Lattice &, const unsigned, const unsigned) const;
  inline void collide_and_bound(Lattice &, IncompFlowMultiscaleMap &,
                                const IncompFlowCollisionManager &,
                                const unsigned, const unsigned) const;
  virtual ~AbstractNodeDesc() = 0;

private:
  virtual void stream_(Lattice &, const unsigned, const unsigned) const
      noexcept = 0;
  virtual void stream_with_bcheck_(Lattice &, const unsigned,
                                   const unsigned) const = 0;
  virtual void collide_and_bound_(Lattice &, IncompFlowMultiscaleMap &,
                                  const IncompFlowCollisionManager &,
                                  const unsigned, const unsigned) const;
};

//! D2Q9 base class collide and bound
//!
//! \param lat Lattice
//! \param mmap Incompressible flow multiscale map
//! \param cman Collision manager
//! \param i Index of node in the y-direction
//! \param j Index of node in the x-direction
inline void
AbstractNodeDesc::collide_and_bound_(Lattice &lat,
                                     IncompFlowMultiscaleMap &mmap,
                                     const IncompFlowCollisionManager &cman,
                                     const unsigned i, const unsigned j) const {
  collide_and_bound_(lat, mmap, cman, i, j);
}

//! \class NodeInactive
//!
//! \brief Inactive node //!
//! Represents an inactive node in the domain
class NodeInactive : public AbstractNodeDesc {
public:
  ~NodeInactive() {}

private:
  virtual void stream_(Lattice &, const unsigned, const unsigned) const
      noexcept {}
  virtual void stream_with_bcheck_(Lattice &, const unsigned,
                                   const unsigned) const {}
  void collide_and_bound_(Lattice &, IncompFlowMultiscaleMap &,
                          const IncompFlowCollisionManager &, const unsigned,
                          const unsigned) const {}
};

//! \class AbstractNodeActive
//!
//! \brief Active node
//!
//! Represents an active node in the domain where streaming and collision
//! occur
class AbstractNodeActive : public AbstractNodeDesc {
public:
  virtual ~AbstractNodeActive() = 0;

protected:
  virtual void stream_(Lattice &, const unsigned, const unsigned) const
      noexcept;
  virtual void stream_with_bcheck_(Lattice &, const unsigned,
                                   const unsigned) const;
  virtual void collide_and_bound_(Lattice &, IncompFlowMultiscaleMap &,
                                  const IncompFlowCollisionManager &,
                                  const unsigned, const unsigned) const;
};

//! \class NodeActive
//!
//! \brief Concrete active node
//!
//! Concrete class of an active node
class NodeActive : public AbstractNodeActive {
public:
  ~NodeActive() {}
};

//! \class NodeWestFacingWall
//!
//! \brief West facing wall
//!
//! Represents a solid, west facing wall
class NodeWestFacingWall : public AbstractNodeDesc {
public:
  ~NodeWestFacingWall() {}

private:
  void stream_(Lattice &, const unsigned, const unsigned) const noexcept;
  void stream_with_bcheck_(Lattice &, const unsigned, const unsigned) const;
  void collide_and_bound_(Lattice &, IncompFlowMultiscaleMap &,
                          const IncompFlowCollisionManager &, const unsigned,
                          const unsigned) const;
};

//! \class NodeSouthFacingWall
//!
//! \brief South facing wall
//!
//! Represents a solid, south facing wall
class NodeSouthFacingWall : public AbstractNodeDesc {
public:
  ~NodeSouthFacingWall() {}

private:
  void stream_(Lattice &, const unsigned, const unsigned) const noexcept;
  void stream_with_bcheck_(Lattice &, const unsigned, const unsigned) const;
  void collide_and_bound_(Lattice &, IncompFlowMultiscaleMap &,
                          const IncompFlowCollisionManager &, const unsigned,
                          const unsigned) const;
};

//! \class NodeEastFacingWall
//!
//! \brief East facing wall
//!
//! Represents a solid, east facing wall
class NodeEastFacingWall : public AbstractNodeDesc {
public:
  ~NodeEastFacingWall() {}

private:
  void stream_(Lattice &, const unsigned, const unsigned) const noexcept;
  void stream_with_bcheck_(Lattice &, const unsigned, const unsigned) const;
  void collide_and_bound_(Lattice &, IncompFlowMultiscaleMap &,
                          const IncompFlowCollisionManager &, const unsigned,
                          const unsigned) const;
};

//! \class NodeNorthFacingWall
//!
//! \brief North facing wall
//!
//! Represents a solid, North facing wall
class NodeNorthFacingWall : public AbstractNodeDesc {
public:
  ~NodeNorthFacingWall() {}

private:
  void stream_(Lattice &, const unsigned, const unsigned) const noexcept;
  void stream_with_bcheck_(Lattice &, const unsigned, const unsigned) const;
  void collide_and_bound_(Lattice &, IncompFlowMultiscaleMap &,
                          const IncompFlowCollisionManager &, const unsigned,
                          const unsigned) const;
};

//! \class Periodic boundary condition
//!
//! \brief Implements a periodic boundary condition
class NodePeriodic : public AbstractNodeActive {
public:
  ~NodePeriodic() {}
  NodePeriodic(const unsigned, const unsigned, const unsigned *,
               const unsigned);

private:
  void collide_and_bound_(Lattice &, IncompFlowMultiscaleMap &,
                          const IncompFlowCollisionManager &, const unsigned,
                          const unsigned) const;
  unsigned i_next_;
  unsigned j_next_;
  std::unique_ptr<unsigned[]> ks_;
  unsigned nk_;
};

//! Constant expression for maximum node descriptor size
//!
//! \return Maximum node descriptor size
constexpr std::size_t max_node_desc_size() {
  return std::max({sizeof(NodeActive), sizeof(NodeWestFacingWall),
                   sizeof(NodeSouthFacingWall), sizeof(NodeEastFacingWall),
                   sizeof(NodeNorthFacingWall)});
}

} // namespace d2q9

} // namespace balbm

#endif // NODE_DESC_HH
