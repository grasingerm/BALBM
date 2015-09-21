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

#include "simulate.hh"

namespace balbm {

namespace d2q9 {

//! Virtual constructor definition
AbstractSimulation::~AbstractSimulation() {}

//! Constructor for incompressible flow simulation
//!
//! \param ni Number of nodes in the y-direction
//! \param nj Number of nodes in the x-direction
//! \param rho Reference density
//! \param mu Reference kinematic viscosity
//! \param pfeq Pointer to base class for incomp flow equilibirum functions
//! \param pconstiteq Pointer to base class for constitutive equations
//! \param pforce Pointer to base class for external forcing scheme
//! \param scbs Vector of callback functions to execute after each time step
IncompFlowSimulation::IncompFlowSimulation(
    const unsigned ni, const unsigned nj, const double rho, const double mu,
    AbstractIncompFlowEqFunct *pfeq, AbstractConstitutiveEq *pconstiteq,
    AbstractForce *pforce, std::vector<AbstractSimCallback *> *pscbs = nullptr)
    : AbstractSimulation(), lat_(Lattice(ni, nj, rho)),
      mmap_(IncompFlowMultiscaleMap(ni, nj, 
            mu_to_omega(mu, lat_.cssq(), lat_.dt()))),
      cman_(IncompFlowCollisionManager(pfeq, pconstiteq, pforce)),
      spscbs_(std::unique_ptr<std::vector<AbstractSimCallback *>>(pscbs)) {}

//! Run an imcompressible flow simulation
//!
//! \param nsteps Steps to simulate
//! \return number of steps simulated
unsigned IncompFlowSimulation::simulate_(const unsigned nsteps) {
  unsigned init_step = step();

  try {
    for (unsigned k = init_step; k <= nsteps; ++k)
      simulate_();
  } catch (std::exception &e) {
    std::cerr << "ERROR: simulation terminated after " << step() << " steps.\n"
              << e.what() << '\n';
    throw;
  }

  return nsteps - init_step;
}

//! Simulate a time step
unsigned IncompFlowSimulation::simulate_() {
  // code for one time step
  lat_.stream();
  lat_.swap_f_ptrs();
  lat_.collide_and_bound(mmap_, cman_);
  if (spscbs_)
    for (const auto &cb : *spscbs_)
      (*cb)(*this);
  ++step_;
  return 1;
}

} // namespace d2q9

} // namespace balbm
