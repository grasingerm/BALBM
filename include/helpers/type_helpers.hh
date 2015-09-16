#ifndef TYPE_HELPERS_HH
#define TYPE_HELPERS_HH

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

//! compiler dianostic for checking a deduced type
//! usage: TypeDeduction<decltype(var)> vartype;
template <typename T> class TypeDeduction;

//! template metaprogramming for removing const from a type
template <class T> struct remove_const { using type = T; };

template <class U> struct remove_const<U const> { using type = U; };

#endif // TYPE_HELPERS_HH
