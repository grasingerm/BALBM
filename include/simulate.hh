#ifndef SIMULATE_HH
#define SIMULATE_HH

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

#include "callback.hh"
#include "collision_manager.hh"
#include "lattice.hh"
#include "multiscale_map.hh"
#include <memory>
#include <vector>

namespace balbm {

namespace d2q9 {

//! \class AbstractSimulation
//!
//! \brief Abstract base class for simulation types
class AbstractSimulation {
public:
  virtual ~AbstractSimulation() = 0;
  AbstractSimulation() : step_(0) {}
  inline unsigned simulate(const unsigned nsteps) { return simulate_(nsteps); }
  inline unsigned step() const { return step_; }

protected:
  unsigned step_;

private:
  virtual unsigned simulate_(const unsigned nsteps) = 0;
};

// TODO: should we push inheritance hierarchies and customization down the
//       composition chain? or up? down is more modular and complex, up is
//       more efficient and rigid... both is the most flexible?

//! \class IncompFlowSimulation
//!
//! \brief Incompressible flow lattice Boltzmann method simulation
class IncompFlowSimulation : public AbstractSimulation {
public:
  ~IncompFlowSimulation() {}
  IncompFlowSimulation(const unsigned, const unsigned, const double,
                       const double, AbstractIncompFlowEqFunct *,
                       AbstractConstitutiveEq *, AbstractForce *,
                       std::vector<AbstractSimCallback *> *);
  inline const IncompFlowMultiscaleMap &multiscale_map() const { return mmap_; }
  template <typename Node, typename... Args>
  inline void set_node_desc(unsigned i, unsigned j, Args... args) {
    lat_.set_node_desc<Node>(i, j, args...);
  }

private:
  unsigned simulate_(const unsigned);
  unsigned simulate_();
  Lattice lat_;
  IncompFlowMultiscaleMap mmap_;
  IncompFlowCollisionManager cman_;
  std::unique_ptr<std::vector<AbstractSimCallback *>> spscbs_;
};

} // namespace d2q9

} // namespace balbm

#endif // SIMULATE_HH
