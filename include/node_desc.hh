#ifndef __NODE_DESC_HH__
#define __NODE_DESC_HH__

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

#include "balbm_config.hh"

namespace balbm
{

namespace d2q9
{
//! \class AbstractNodeDesc
//!
//! \brief Abstract base class for  node descriptors
//!
//! Provide polymorphic behavior for each node in the lattice based on its
//! physical "status" as a node, e.g. change behavior of streaming and collision
//! steps in order to simulate appropriate physics and boundary conditions
class AbstractNodeDesc
{
public:
  inline void stream(Lattice&, const unsigned, const unsigned) const;
  inline void collide_and_bound(Lattice&, const CollisionManager&,
                                const unsigned, const unsigned) const noexcept;
  virtual ~NodeDesc()=0;
private:
  virtual void _stream(Lattice&, const unsigned, const unsigned) const
                       noexcept=0;
  virtual void _stream_with_bcheck(Lattice&, const unsigned, 
                                   const unsigned) const=0;
  virtual void _collide_and_bound(Lattice&, const CollisionManager&,
                                  const unsigned, const unsigned) 
                                  const noexcept=0;
};

//! \class NodeInactive
//!
//! \brief Inactive node //!
//! Represents an inactive node in the domain
class NodeInactive : public AbstractNodeDesc
{
public:
  ~NodeInactive() {}
private:
  virtual void _stream(Lattice&, const unsigned, const unsigned) const
                       noexcept {}
  virtual void _stream_with_bcheck(Lattice&, const unsigned, 
                                   const unsigned) const {}
  virtual void _collide_and_bound(Lattice&, const CollisionManager&,
                                  const unsigned, const unsigned) 
                                  const noexcept {};
};

//! \class AbstractNodeActive
//!
//! \brief Active node
//!
//! Represents an active node in the domain where streaming and collision
//! occur
class AbstractNodeActive : public AbstractNodeDesc
{
public:
  virtual ~NodeActive()=0;
private:
  virtual void _stream(Lattice&, const unsigned, const unsigned) const noexcept;
  virtual void _stream_with_bcheck(Lattice&, const unsigned, 
                           const unsigned) const;
  virtual void _collide_and_bound(Lattice&, const CollisionManager&,
                          const unsigned, const unsigned) const noexcept;
};

//! \class NodeActive
//!
//! \brief Concrete active node
//!
//! Concrete class of an active node
class NodeActive : public AbstractNodeActive
{
public:
  ~NodeActive() {}
};

//! \class NodeWestFacingWall
//!
//! \brief West facing wall
//!
//! Represents a solid, west facing wall
class NodeWestFacingWall : public AbstractNodeActive
{
public:
  ~NodeWestFacingWall() {}
private:
  void _stream(Lattice&, const unsigned, const unsigned) const;
  void _stream_with_bcheck(Lattice&, const unsigned, 
                           const unsigned) const;
  void _collide_and_bound(Lattice&, const CollisionManager&,
                          const unsigned, const unsigned) const;
};

//! \class NodeSouthFacingWall
//!
//! \brief South facing wall
//!
//! Represents a solid, south facing wall
class NodeSouthFacingWall : public AbstractNodeActive
{
public:
  ~NodeSouthFacingWall() {}
private:
  void _stream(Lattice&, const unsigned, const unsigned) const;
  void _stream_with_bcheck(Lattice&, const unsigned, 
                           const unsigned) const;
  void _collide_and_bound(Lattice&, const CollisionManager&,
                          const unsigned, const unsigned) const;
};

//! \class NodeEastFacingWall
//!
//! \brief East facing wall
//!
//! Represents a solid, east facing wall
class NodeEastFacingWall : public AbstractNodeActive
{
public:
  ~NodeEastFacingWall() {}
private:
  void _stream(Lattice&, const unsigned, const unsigned) const;
  void _stream_with_bcheck(Lattice&, const unsigned, 
                           const unsigned) const;
  void _collide_and_bound(Lattice&, const CollisionManager&,
                          const unsigned, const unsigned) const;
};

//! \class NodeNorthFacingWall
//!
//! \brief North facing wall
//!
//! Represents a solid, North facing wall
class NodeNorthFacingWall : public AbstractNodeActive
{
public:
  ~NodeNorthFacingWall() {}
private:
  void _stream(Lattice&, const unsigned, const unsigned) const;
  void _stream_with_bcheck(Lattice&, const unsigned, 
                           const unsigned) const;
  void _collide_and_bound(Lattice&, const CollisionManager&,
                          const unsigned, const unsigned) const;
};

} // namespace d2q9

} // namespace balbm

#endif // __NODE_DESC_HH__
