cmake_minimum_required(VERSION 3.14)
project(orbital_propagator)

set(CMAKE_CXX_STANDARD 14)

set(PROJECT_SOURCE_DIR "/Users/sylvain/Developer/orbital_propagator")

include_directories(. ../libs ../libs/boost_1_69_0 ../libs/cspice/include ../libs/earthGRAM2016/src)
#link_directories(../libs/earthGRAM2016/lib ../libs/cspice/lib)

#add_executable(orbital_propagation

        #)

#target_link_libraries(orbital_propagation libearthGRAM2016.a cspice.a csupport.a)

add_library(libearthGRAM2016 STATIC IMPORTED) # or STATIC instead of SHARED
set_target_properties(libearthGRAM2016 PROPERTIES
        IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/libs/earthGRAM2016/lib/libearthGRAM2016.a"
        INTERFACE_INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/libs/earthGRAM2016/src"
        )

add_library(cspice STATIC IMPORTED) # or STATIC instead of SHARED
set_target_properties(cspice PROPERTIES
        IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/libs/cspice/lib/cspice.a"
        INTERFACE_INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/libs/cspice/include"
        )

add_library(csupport STATIC IMPORTED) # or STATIC instead of SHARED
set_target_properties(csupport PROPERTIES
        IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/libs/cspice/lib/csupport.a"
        INTERFACE_INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/libs/cspice/include"
        )

set(ORBITAL_PROPAGATOR_SRCS
        constants.h
        EnvironmentModel.cpp
        EnvironmentModel.h
        EqMotion.cpp
        EqMotion.h
        main.cpp
        perturbations.cpp
        perturbations.h
        Spacecraft.cpp
        Spacecraft.h
        state_type.h
        )

add_executable(orbital_propagator ${ORBITAL_PROPAGATOR_SRCS})
target_link_libraries(orbital_propagator libearthGRAM2016 cspice csupport)