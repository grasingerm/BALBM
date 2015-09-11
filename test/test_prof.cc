#include "helpers/prof_helpers.hh"
#include <thread>

using namespace baprof;

int main() {
  tic();
  std::this_thread::sleep_for(std::chrono::milliseconds(3000));
  toc();

  std::cout << "Should pause and display approximately 3 seconds\n";

  return 0;
}
