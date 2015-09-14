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

#ifndef BALBM_CONFIG_HH
#define BALBM_CONFIG_HH

// NOTE: for maximum performance define NDEBUG and do not define
//       BALBM_CHECK_BOUNDS_STREAMING.

//#define NDEBUG

//! Define this to use runtime bounds checking when streaming
//#define BALBM_CHECK_BOUNDS_STREAMING

//! Define this to disable armadillo bounds checking
//#define ARMA_NO_DEBUG

//! Define this to allow armadillo to save to hdf5
//#define ARMA_USE_HDF5

namespace balbm {
const static char *VERSION = "0.0.1";
}

#endif // BALBM_CONFIG_HH
