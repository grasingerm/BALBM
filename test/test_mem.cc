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

#include "balbm.hh"
#include <cassert>
#include <iostream>

using namespace balbm::d2q9;
using namespace std;

//! simulation parameters
const static unsigned ni = 40;
const static unsigned nj = 12;

int main() {
  cout << "Maximum node descriptor size = " << max_node_desc_size() << '\n';
  SimpleMemPool mem_pool(ni * nj * max_node_desc_size());
  cout << "Memory capacity = " << mem_pool.capacity() << '\n';
  cout << "What should the memory capacity be? = "
       << ni *nj *max_node_desc_size() << '\n';

  return 0;
}
