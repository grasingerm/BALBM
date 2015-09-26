cmake_minimum_required(VERSION 2.8)

project(BALBM)

# dependencies
add_executable(test_mem test_mem.cc)
add_executable(test_prof test_prof.cc)
add_executable(test_lat_vecs test_lat_vecs.cc ../src/lattice.cc)
add_executable(test_poiseuille_newtonian test_poiseuille_newtonian.cc
                                         ../src/collision_manager.cc
                                         ../src/constitutive.cc
                                         ../src/equilibrium.cc
                                         ../src/force.cc
                                         ../src/lattice.cc
                                         ../src/multiscale_map.cc
                                         ../src/node_desc.cc
                                         ../src/simulate.cc             )

# link libraries
target_link_libraries(test_poiseuille_newtonian armadillo)
if (${UNIX})
  target_link_libraries(test_poiseuille_newtonian m)
endif ()


# install
install(
        TARGETS 
                test_mem
                test_prof
                test_lat_vecs
                test_poiseuille_newtonian
        DESTINATION 
                tests
       )
