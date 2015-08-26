#ifndef __NODE_DESC_HH__
#define __NODE_DESC_HH__

//! \class NodeDescD2Q9
//!
//! \brief Abstract base class for D2Q9 node descriptors
//!
//! Provide polymorphic behavior for each node in the lattice based on its
//! physical "status" as a node, e.g. change behavior of streaming and collision
//! steps in order to simulate appropriate physics and boundary conditions
class NodeDescD2Q9
{
public:
  virtual void stream(LatticeD2Q9&, unsigned i, unsigned j)=0;
  virtual ~NodeDescD2Q9();
};



#endif // __NODE_DESC_HH__
