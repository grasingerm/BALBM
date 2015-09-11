#ifndef MULTISCALE_MAP_HH
#define MULTISCALE_MAP_HH

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

#include "lattice.hh"
#include <memory>

namespace balbm {

namespace d2q9 {

//! Convert visocisty to relaxation time
//!
//! \param mu Kinematic viscosity
//! \param cssq Speed of sound squared
//! \param dt Length of time step
//! \return Relaxation time
inline double mu_to_relax(const double mu, const double cssq, const double dt) {
  return mu / (cssq * dt) + 0.5;
}

//! Convert visocisty to collision frequency
//!
//! \param mu Kinematic viscosity
//! \param cssq Speed of sound squared
//! \param dt Length of time step
//! \return Collision frequency
inline double mu_to_omega(const double mu, const double cssq, const double dt) {
  return 1.0 / mu_to_relax(mu, cssq, dt);
}

//! \class AbstractMultiscaleMap
//!
//! \brief Base class for map from mesoscale to macroscale
//!
//! Maps particle distributions to macroscopic variables of interest
class AbstractMultiscaleMap {
public:
  AbstractMultiscaleMap(const unsigned nx, const unsigned ny)
      : ni_(ni), nj_(nj), rho_(std::unique_ptr<double[]>(new double[nx * ny])) {
  }
  virtual ~AbstractMultiscaleMap() = 0;
  inline double num_i() const noexcept { return ni_; }
  inline double num_j() const noexcept { return nj_; }
  inline double rho(const unsigned i, const unsigned j) const {
    return sprho_[i * num_j() + j];
  }
  inline void map_to_macro(const Lattice &lat) { map_to_macro_(lat); }

protected:
  inline double &rho_(const unsigned i, const unsigned j) {
    return sprho_[i * num_j() + j];
  }
  virtual void map_to_macro_(const Lattice &);
  virtual void map_to_macro_(const Lattice &, const unsigned, const unsigned);

private:
  unsigned ni_;
  unsigned nj_;
  std::unique_ptr<double[]> sprho_;
};

//! \class DensityMultiscaleMap
//!
//! \brief Maps particle distributions to local densities
//!
//! Concrete class for density-based multiscale map
class DensityMultiscaleMap : public AbstractMultiscaleMap {
public:
  DensityMultiscaleMap(const unsigned ni, const unsigned nj)
      : AbstractMultiscaleMap(ni, nj) {}
  ~DensityMultiscaleMap() {}
};

//! \class IncompFlowMultiscaleMap
//!
//! \brief Maps particle distributions to local macroscopic flow variables
//!
//! Concrete class for incompressible flow multiscale map. Maps particle
//! distributions to local macroscopic density, flow, and collision frequency
class IncompFlowMultiscaleMap : public AbstractMultiscaleMap {
public:
  IncompFlowMultiscaleMap(const unsigned ni, const unsigned nj,
                          const double omega)
      : AbstractMultiscaleMap(ni, nj),
        spu_(std::unique_ptr<double[]>(new double[ni * nj * 2])),
        spomega_(std::unique_ptr<double[]>(new double[ni * nj])) {
    init_(omega)
  }
  ~IncompFlowMultiscaleMap() {}
  inline double u(const unsigned i, const unsigned j, const unsigned c) const {
    return spu_[2 * (i * num_j() + j) + c];
  }
  inline const double *pu(const unsigned i, const unsigned j) const {
    return &spu_[2 * (i * num_j() + j)];
  }
  inline double &omega(const unsigned i, const unsigned j) {
    return spomega_[i * num_j() + j];
  }
  inline double omega(const unsigned i, const unsigned j) const {
    return spomega_[i * num_j() + j];
  }

private:
  inline double &u_(const unsigned i, const unsigned j, const unsigned c) {
    return spu_[2 * (i * num_j() + j) + c];
  }
  void map_to_macro_(const Lattice &, const unsigned, const unsigned);
  void init_();
  std::unique_ptr<double[]> spu_;
  std::unique_ptr<double[]> spomega_;
};

} // namespace d2q9

} // namespace balbm

#endif // MULTISCALE_MAP_HH
