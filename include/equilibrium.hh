#ifndef EQUILIBRIUM_HH
#define EQUILIBRIUM_HH

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

#include <armadillo>

namespace balbm {

namespace d2q9 {

/* TODO: does it make sense to have a base class for all equilibrium equations?
//! \class AbstractEqFunct
//!
//! \brief Abstract class for equilibrium distribution function
//!
//! Defines functor for calculating a local equilibrium distribution function
class AbstractEqFunct
{
public:
  virtual ~AbstractEqFunct()=0;
  inline double f const(const Lattice& lat, const AbstractMultiscaleMap& mmap,
                        const unsigned k)
    { return f_(lat, mmap, k); }
private:
  virtual double f_ const(const Lattice& lat, const AbstractMultiscaleMap& mmap,
                          const unsigned k)=0;
};
*/

// TODO: if we start implementing other lattices then "f" should be a function
// template/overloaded?

//! \class AbstractIncompFlowEqFunct
//!
//! \brief Base equilibrium distribution function for incompressible flow
//!
//! Defines functor for calculating a local equilibrium distribution function
//! for simulating incompressible flow
class AbstractIncompFlowEqFunct {
public:
  virtual ~AbstractIncompFlowEqFunct() = 0;
  inline double f(const Lattice &lat, const double rho, const arma::vec &u,
                  const unsigned k) const {
    return f_(lat, rho, u, k);
  }

private:
  virtual double f_(const Lattice &lat, const double rho, const arma::vec &u,
                    const unsigned k) const = 0;
};

//! \class IncompFlowEqFunct
//!
//! \brief Equilibrium distribution function for incompressible flow
//!
//! Defines functor for calculating a local equilibrium distribution function
//! for simulating incompressible flow
class IncompFlowEqFunct : public AbstractIncompFlowEqFunct {
public:
  ~IncompFlowEqFunct() {}

private:
  double f_(const Lattice &lat, const double rho, const arma::vec &u,
            const unsigned k) const;
};

//! \class IncompFlowHLEqFunct
//!
//! \brief He and Lou equilibrium distribution function for incompressible flow
//!
//! Defines functor for calculating a local He and Lou equilibrium distribution
//! function for simulating incompressible flow. Used for numerical stability.
class IncompFlowHLEqFunct : public AbstractIncompFlowEqFunct {
public:
  ~IncompFlowHLEqFunct() {}
  IncompFlowHLEqFunct(unsigned rho_o) : rho_o_(rho_o) {}

private:
  double rho_o_;
  double f_(const Lattice &lat, const double rho, const arma::vec &u,
            const unsigned k) const = 0;
};

} // namespace d2q9

} // namespace balbm

#endif // EQUILIBRIUM_HH
