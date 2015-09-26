#include "lattice.hh"
#include <cassert>
#include <iostream>

using namespace balbm::d2q9;
using namespace std;

int main()
{
  Lattice lat;

  cout << "Testing public interface for lattice vectors...\n";

  assert(1 == lat.c(1, 0));
  assert(0 == lat.c(1, 1));
  assert(0 == lat.c(2, 0));
  assert(1 == lat.c(2, 1));
  assert(-1 == lat.c(3, 0));
  assert(0 == lat.c(3, 1));
  assert(0 == lat.c(4, 0));
  assert(-1 == lat.c(4, 1));
  assert(1 == lat.c(5, 0));
  assert(1 == lat.c(5, 1));
  assert(-1 == lat.c(6, 0));
  assert(1 == lat.c(6, 1));
  assert(-1 == lat.c(7, 0));
  assert(-1 == lat.c(7, 1));
  assert(1 == lat.c(8, 0));
  assert(-1 == lat.c(8, 1));

  cout << "TEST PASSED\n";

  return 0;
}
