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
#include <cmath>
#include <vector>

using namespace balbm::d2q9;
using namespace std;

//! simulation parameters
const static unsigned               ni              = 40;
const static unsigned               nj              = 12;
const static double                 rho             = 1.0;
const static double                 mu              = 1.0/6.0;
const static double                 pgrad           = -1.102e-3;
const static std::array<double, 2>  F               { { -pgrad, 0.0 } };
const static unsigned               nsteps          = 5000; 
const static unsigned               stepout         = 50; 

int main()
{
  // ... new SukopThorneForce(F)
  // ... new NewtonianConstitutiveEq(mu)
  // ... new IncompFlowEqFunct()
  // ... ... new IncompFlowCollisionManager(new IncompFlowEqFunct(), 
  //                                        new NewtonianConstitutiveEq(mu),
  //                                        new SukopThorneForce(F));
  // ... new Lattice(ni, nj, rho)
  // ... new IncompFlowMultiscaleMap(ni, nj, mu_to_omega(mu))
 
  // TODO: add this to a list of test helper functions
  const auto analytic_soln = [&](const vector<double>& xs) {
    const double h = (nj - 1) / 2.0;
    vector<double> result(xs.size());
    for (unsigned i = 0; i < xs.size(); ++i)
      result[i] = -1.0 / (2.0 * mu) * pgrad * (h*h - x[i]*x[i]);
    return result;
  }

  SimpleMemPool mem_pool;
  mem_pool.reserve(max_sim_callback_size());

  vector<AbstractSimCallback*> scbs;
  scbs.push_back(mem_pool.allocate(DisplayTimestepCallback(stepout));

  // Simulation(Lattice*, IncompFlowMultiscaleMap*, IncompFlowCollisionManager,
  //            vector<AbstractSimCallbacks*>)
  Simulation sim
    (ni, nj, rho, mu_to_omega(mu)),
     new IncompFlowEqFunct(), 
     new NewtonianConstitutiveEq(mu),
     new SukopThorneForce(F),
     &scbs
    );

  baprof::tic();
  sim.simulate(nsteps); // run simulation
  baprof::toc(); // time simulation

  const unsigned i = ni / 2; 
  const auto& mmap = sim.multiscale_map();
  vector<double> xs(nj);
  for (unsigned j = 0; j < nj; ++j) xs[j] = (j - nj/2.0 - 0.5);
  const auto& us = analytic_soln(xs);
  
  for (unsigned j = 0; j < nj; ++j)
  {
    cout << "analyt == lbm ? " << us(j) << " == " << mmap.u(i, j, 0);
    assert(fabs(us(j) - mmap.u(i, j, 0)) / us(j) <= 5e-3);
  }

  cout << "TEST PASSED\n";

  // manually call destructors for mem pool
  // TODO: this isn't exception safe... I don't think
  // TODO: find a better mem pool
  for (auto& scb : scbs) scb->~AbstractSimCallback();

  return 0;
}
