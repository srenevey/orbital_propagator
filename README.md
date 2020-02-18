This orbital propagator computes the trajectory of a spacecraft in the two-body problem with perturbations from atmospheric drag (exponentional model or EarthGRAM), geopotential perturbations, third-body effects, and solar radiation pressure. An example is provided in main.cpp and the documentation can be built by running `doxygen` from `doc`.

Note: the propagator has not been validated. The output might be inaccurate / incorrect.



## Downloading the required models

The C++ propagator uses the EGM2008 geopotential model, EarthGRAM2016 for the atmospheric model as well as several kernels used with the Naif SPICE toolbox for the planetary ephemerides.

### Geopotential model
The EGM2008 model can be downloaded [here](http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/first_release.html). Download the model *Spherical Harmonic Coefficients for Earth's Gravitational Potential - "Tide Free" system* and then extract it. Place the file *EGM2008_to2190_TideFree* into *assets*.

```
orbital_propagator
└── assets
    |   EGM2008_to2190_TideFree
|   doc
|   libs
|   README.md
|   src
```

### Atmosphere model
EarthGRAM2016 can be requested from the [NASA Software Catalog](https://software.nasa.gov/software/MFS-32780-2). Unzip the archive and place it in the *libs* folder

```
orbital_propagator
└── assets
    |   EGM2008_to2190_TideFree 	
|   doc
└── libs
    |   EarthGRAM2016
|   README.md
|   src
```

Copy the file NameRef_Linux.txt from IOfiles into EarthGRAM2016 and rename it to NameRef.txt. Open NameRef.txt and modifiy as follows:

    atmpath = <ABSOLUTE_PATH_TO_PROJECT_ROOT>/libs/EarthGRAM2016/IOfiles/
    NCEPpath = <ABSOLUTE_PATH_TO_PROJECT_ROOT>/libs/EarthGRAM2016/NCEPdata/FixedBin/
    trapath = null
    prtpath = null
    nprpath =  null
    conpath =  null
    rrapath = <ABSOLUTE_PATH_TO_PROJECT_ROOT>/libs/EarthGRAM2016/RRAdata/
    ...
    iurra = 0
    ...
    ibltest = 0

where <ABSOLUTE_PATH_TO_PROJECT_ROOT> must be modified accordingly. By default, EarthGRAM2016 prompts the user to enter the path to NameRef.txt which is not ideal for a seamless integration. The file *earthgram.patch* modifies some of the source files to suppress any prompt and output. To apply the patch, copy the file *earthgram.patch* from *doc* into *libs/EarthGRAM2016/Source* and then in a terminal, navigate to *libs/EarthGRAM2016/Source* and type

    patch -p1 < earthgram.patch

The EarthGRAM2016 model needs to be compiled into a library. Assuming that [CMake](https://cmake.org/) is installed, create a file called CMakeList.txt in EarthGRAM2016/Source and add the following lines

    cmake_minimum_required(VERSION 3.12)
    project(EarthGRAM2016)
    
    set(CMAKE_CXX_STANDARD 11)
    
    add_library(EarthGRAM2016 STATIC
        Atmod1.cpp
        Atmod1.h
        AuxProf.cpp
        AuxProf.h
        HWM.cpp
        HWM.h
        Init.cpp
        Init.h
        InitP.cpp
        InitP.h
        JB2008.cpp
        JB2008.h
        Map.cpp
        Map.h
        MET.cpp
        MET.h
        MSIS.cpp
        MSIS.h
        NCEP.cpp
        NCEP.h
        Pert.cpp
        Pert.h
        RRA.cpp
        RRA.h)
    
    set_property(TARGET EarthGRAM2016 PROPERTY POSITION_INDEPENDENT_CODE ON)

Open a terminal and navigate to EarthGRAM2016/Source. Then, type

    mkdir build
    cd build
    cmake ..
    make

A file EarthGRAM2016.a should be crated in the *build* directory.

## Planetary ephemerides
The high precision propagator uses JPL's ephemerides to retrieve the planet's positions. This process is facilitated by the SPICE toolbox developed by NASA's Naif. The kernels used by the SPICE toolbox can be accessed from the [JPL's Naif website](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/) by going to Data > Generic Kernels > Generic Kernels and navigating in the subfolders. The required kernels are (right click > "Save Link As..." to download the file directly):
- [de430.bsp](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp)
- [earth_latest_high_prec.bpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc)
- [gm_de431.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc)
- [naif0012.tls](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls)
- [pck00010.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc)

These kernels should be placed in a *kernels* folder inside the *assets* directory.

```
orbital_propagator
└── assets
    |   EGM2008_to2190_TideFree 
└── kernels
    |   de430.bsp
    |   earth_latest_high_prec.bpc
    |   gm_de431.tpc
    |   naif0012.tls
    |   pck00010.tpc
|   doc
└── libs
    |   EarthGRAM2016
|   README.md
|   src
```

These kernels will be loaded by SPICE through a meta-kernel. The following content should be put in a file named *kernels.tm* and placed in the *assets* folder. Modify <ABSOLUTE_PATH_TO_PROJECT_ROOT> to reflect your project's root directory path:

    KPL/MK
       \begindata
       PATH_VALUES     = ( '<ABSOLUTE_PATH_TO_PROJECT_ROOT>' )
       
       PATH_SYMBOLS    = ( 'ROOT' )
    
       KERNELS_TO_LOAD = (  '$ROOT/assets/kernels/de430.bsp',
                            '$ROOT/assets/kernels/earth_latest_high_prec.bpc',
                            '$ROOT/assets/kernels/gm_de431.tpc',
                            '$ROOT/assets/kernels/naif0012.tls',
                            '$ROOT/assets/kernels/pck00010.tpc')
     
       \begintext



## Compilation

The orbital propagator is based on Eigen, boost and CSPICE. These libraries should be downloaded and installed before compiling the orbital propagator. An example of CMakeLists.txt file is shown below, assuming that the libraries are located in *orbital_propagator/libs*:

```cmake
cmake_minimum_required(VERSION 3.14)
project(orbital_propagator)

set(CMAKE_CXX_STANDARD 14)

set(PROJECT_SOURCE_DIR "<ORBITAL_PROPAGATOR_ROOT_DIR>")

include_directories(. ../libs ../libs/boost_1_69_0 ../libs/cspice/include ../libs/earthGRAM2016/src)

add_library(libearthGRAM2016 STATIC IMPORTED)
set_target_properties(libearthGRAM2016 PROPERTIES
        IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/libs/earthGRAM2016/lib/libearthGRAM2016.a"
        INTERFACE_INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/libs/earthGRAM2016/src"
        )

add_library(cspice STATIC IMPORTED)
set_target_properties(cspice PROPERTIES
        IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/libs/cspice/lib/cspice.a"
        INTERFACE_INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/libs/cspice/include"
        )

add_library(csupport STATIC IMPORTED)
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
```



## Building the doc

Install doxygen and then from the *doc* directory enter

```bash
doxygen
```

The main entry point is *doc/html/index.html*.

