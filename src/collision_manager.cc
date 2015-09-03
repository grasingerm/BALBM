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
#include <armadillo>

namespace balbm
{

namespace d2q9
{

//! Virtual destructor
AbstractCollisionManager::~AbstractCollisionManager() {}

//! Incompressible flow collision
//!
//! \param lat Lattice
//! \param mmap Multiscale map
//! \param i Index in x-direction
//! \param j Index in y-direction
void IncompFlowCollisionManager::collide_
  (Lattice& lat, IncompFlowMultiscaleMap& mmap, const unsigned i, 
   const unsigned j) const
{
  const auto rhoij = mmap.rho(i, j);
  arma::vec::fixed<2> uij = arma::vec(mmap.pu(i, j));
  if (pextforce_ != nullptr) uij = pextforce_->u_trans(uij);

  const unsigned nk = lat.num_k();
  arma::vec feq(nk);
  arma::vec fneq(nk);
  for (unsigned k = 0; k < nk; ++k)
  {
    feq(k) = pfeq_->f(lat, rhoij, uij, k);
    fneq(k) = lat.f(i, j, k) - feq(k);
  }

  const auto mu = pconstiteq_->mu(lat, mmap, fneq, i, j);
  const auto omega = mu_to_omega(mu, lat.cssq(), lat.dt());
  
  if (pextforce_ != nullptr)
    for (unsigned k = 0; k < nk; ++k)
      lat.f(i, j, k) = (omega * feq(k) + (1.0 - omega) * lat.k(i, j, k)
                        + pextforce_->f_col(lat, omega, uij, k);
  else
    for (unsigned k = 0; k < nk; ++k)
      lat.f(i, j, k) = (omega * feq(k) + (1.0 - omega) * lat.k(i, j, k);

  mmap.omega(i, j) = omega; 
}

} // namespace d2q9

} // namespace balbm
