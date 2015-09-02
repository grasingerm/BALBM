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
#include "lattice.hh"

namespace balbm
{

namespace d2q9
{

//! \class AbstractCollisionManager
//!
//! \brief Abstract base class for collision managers
//!
//! Wrapper class for combining behavior from a constitutive relationship,
//! a solver, and external forces to perform an LBM collision
class AbstractCollisionManager
{
public:
  virtual ~AbstractCollisionManager()=0;
  void collide(Lattice& lat, AbstractMultiscaleMap& mmap) const
    { return collide_(lat, mmp); }
private:
  virtual void collide_(Lattice&, AbstractMultiscaleMap&) const=0;
};

} // namespace d2q9

} // namespace balbm

#endif // COLLISION_MANAGER_HH
