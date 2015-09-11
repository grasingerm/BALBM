#ifndef COLLISION_MANAGER_HH
#define COLLISION_MANAGER_HH

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
#include "constitutive.hh"
#include "lattice.hh"
#include <memory>

namespace balbm {

namespace d2q9 {

// TODO: consider making better use of return value optimization
// TODO: all customizability may come from composition.
//       Do we even need inheritance here?

//! \class AbstractCollisionManager
//!
//! \brief Abstract base class for collision managers
//!
//! Wrapper class for combining behavior from a constitutive relationship,
//! a solver, and external forces to perform an LBM collision
class AbstractCollisionManager {
public:
  virtual ~AbstractCollisionManager() = 0;
  inline void collide(Lattice &lat, AbstractMultiscaleMap &mmap,
                      const unsigned i, const unsigned j) const {
    return collide_(lat, mmp, i, j);
  }

private:
  virtual void collide_(Lattice &, AbstractMultiscaleMap &, unsigned,
                        unsigned) const = 0;
};

//! \class IncompFlowCollisionManager
//!
//! \brief Collision manager for incompressible flow
class IncompFlowCollisionManager {
public:
  ~IncompFlowCollisionManager() {}
  IncompFlowCollisionManager(const AbstractIncompFlowEqFnct *aef,
                             const AbstractConsitutiveEq *ace,
                             const AbstractForce *af = nullptr)
      : pfeq_(std::unique_ptr<AbstractIncompFlowEqFnct>(aef)),
        pconstiteq_(std::unique_ptr<AbstractConstitutiveEq>(ace)),
        pextforce_(std::unique_ptr<AbstractForce>(af)) {}
  void collide_(Lattice &, IncompFlowMultiscaleMap &, const unsigned,
                const unsigned) const;

private:
  std::unique_ptr<AbstractEquilibriumFnct> pfeq_;
  std::unique_ptr<AbstractConstitutiveEq> pconstiteq_;
  std::unique_ptr<AbstractForce> pextforce_;
};

} // namespace d2q9

} // namespace balbm

#endif // COLLISION_MANAGER_HH
