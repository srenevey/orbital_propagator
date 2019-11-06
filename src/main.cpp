#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <boost/numeric/odeint.hpp>
#include "state_type.h"
#include "EnvironmentModel.h"
#include "EqMotion.h"
#include "perturbations.h"

extern "C" {
    #include "SpiceUsr.h"
};

// Use a fifth order Dormand-Prince integration scheme
typedef boost::numeric::odeint::runge_kutta_dopri5< state_type > stepper_type;

// Define a structure to save the time and states
struct push_back_state_and_time
{
    std::vector< state_type >& states_;
    std::vector< double >& times_;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
            : states_( states ) , times_( times ) { }

    void operator()( const state_type &x , double t )
    {
        states_.push_back( x );
        times_.push_back( t );
    }
};


int main() {

    // Define paths to the meta-kernel, geopotential model and EarthGRAM lib
    const char* meta_kernel = "<PORJECT_ROOT>/assets/kernels.tm";
    std::string geopot_model_path = "<PORJECT_ROOT>/assets/EGM2008_to2190_TideFree";
    std::string earthgram_path = "<PORJECT_ROOT>/libs/earthGRAM2016/";


	// Load SPICE kernels
    furnsh_c(meta_kernel);


    // Define the initial epoch, final time (sec), and time interval between to states (sec) for saving purpose, not integration.
    std::string epoch = "2019-03-06 16:00:00";
    double t1 = 48 * 3600.0;
    const double dt = 120.0;

    // Create a spacecraft object
    Spacecraft sc(300.0, 2.1, 3.0, 0.88, 6.0);

    // Create the environment model
    // Third-body flags are: Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune
    std::array<int, 10> third_body_flags = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0};
    EnvironmentModel env_model(earthgram_path, "Earth", 8, geopot_model_path, third_body_flags, true, "EarthGRAM", true, epoch);


    // Define the initial state. Units are km and km/s
    state_type initial_state = {7009.22977, 0.0, 0.0, 0.0, 7.0033, 3.2657};


    SpiceDouble et;
    str2et_c(epoch.c_str(), &et);
    EqMotion system(sc, env_model, et);


    // Integrate the equations of motion.
    std::vector<state_type> states;
    std::vector<double> times;
    double t0 = 0.0;
    boost::numeric::odeint::integrate_const( boost::numeric::odeint::make_controlled( 1e-12 , 1e-12 , stepper_type() ) , system , initial_state, t0, t1, dt, push_back_state_and_time(states, times));


    // Clean up kernel space
    unload_c(meta_kernel);

    // Print the final position
    state_type final_state = states.back();
    std::cout << "Final position (km): \nx: " << std::setprecision(8) << final_state[0] << "\n" << "y: " << std::setprecision(8) << final_state[1] << "\n" << "z: " << std::setprecision(8) << final_state[2] << std::endl;

    return 0;
}