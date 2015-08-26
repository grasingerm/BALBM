#ifndef __NODE_DESC_HH__
#define __NODE_DESC_HH__

namespace balbm
{

namespace d2q9
{
//! \class NodeDesc
//!
//! \brief Abstract base class for  node descriptors
//!
//! Provide polymorphic behavior for each node in the lattice based on its
//! physical "status" as a node, e.g. change behavior of streaming and collision
//! steps in order to simulate appropriate physics and boundary conditions
class NodeDesc
{
public:
  virtual void stream const(Lattice&, const unsigned, const unsigned)=0;
  virtual void collide_and_bound const (Lattice&, const CollisionManager&,
                                        const unsigned, const unsigned)=0;
  virtual ~NodeDesc()=0;
};

//! \class NodeDescInactive
//!
//! \brief Inactive node
//!
//! Represents an inactive node in the domain
class NodeInactive : public NodeDesc
{
public:
  inline void stream const (Lattice& lat, 
                            const unsigned i, const unsigned j) {}
  inline void collide_and_bound const 
    (Lattice& lat, const CollisionManager& cman, const unsigned i,
     const unsigned j) {}
};

//! \class NodeDescActive
//!
//! \brief Active node
//!
//! Represents an active node in the domain where streaming and collision
//! occur
class NodeActive : public NodeDesc
{
public:
  virtual void stream const (Lattice&, const unsigned, const unsigned);
  virtual void collide_and_bound const (Lattice&, const CollisionManager&,
                                        const unsigned, const unsigned);
  virtual ~NodeActive() {}
};

//! \class NodeWestFacingWall
//!
//! \brief West facing wall
//!
//! Represents a solid, west facing wall
class NodeWestFacingWall : public NodeActive
{
public:
  void stream const (Lattice&, const unsigned, const unsigned);
  void collide_and_bound const (Lattice&, const CollisionManager&,
                                        const unsigned, const unsigned);
  ~NodeWestFacingWall() {}
};

//! \class NodeSouthFacingWall
//!
//! \brief South facing wall
//!
//! Represents a solid, south facing wall
class NodeSouthFacingWall : public NodeActive
{
public:
  void stream const (Lattice&, const unsigned, const unsigned);
  void collide_and_bound const (Lattice&, const CollisionManager&,
                                        const unsigned, const unsigned);
  ~NodeSouthFacingWall() {}
};

//! \class NodeEastFacingWall
//!
//! \brief East facing wall
//!
//! Represents a solid, east facing wall
class NodeEastFacingWall : public NodeActive
{
public:
  void stream const (Lattice&, const unsigned, const unsigned);
  void collide_and_bound const (Lattice&, const CollisionManager&,
                                        const unsigned, const unsigned);
  ~NodeEastFacingWall() {}
};

//! \class NodeNorthFacingWall
//!
//! \brief North facing wall
//!
//! Represents a solid, North facing wall
class NodeNorthFacingWall : public NodeActive
{
public:
  void stream const (Lattice&, const unsigned, const unsigned);
  void collide_and_bound const (Lattice&, const CollisionManager&,
                                        const unsigned, const unsigned);
  ~NodeNorthFacingWall() {}
};

} // namespace d2q9

} // namespace balbm

#endif // __NODE_DESC_HH__
