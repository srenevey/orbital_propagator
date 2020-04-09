//
// Created by Sylvain  on 3/31/20.
//

#ifndef ORBITAL_PROPAGATOR_SIM_H
#define ORBITAL_PROPAGATOR_SIM_H


#include "EnvironmentModel.h"
#include "dimensions/Dimensions.h"
#include "Integrator.h"
extern "C" {
    #include "SpiceUsr.h"
};

/** Base class used to load the SPICE kernels and run the simulation. */
class Sim {
    string m_meta_kernel;
    EnvironmentModel m_env_model;

public:
    explicit Sim(std::string meta_kernel);
    ~Sim();

    /** Integrates the equations of motion. Currently uses a Runge-Kutta 4 integration scheme with fixed time step.
     *
     * @param sc            Reference to the spacecraft
     * @param env_model     Environment model
     * @param epoch         Initial epoch with format "YYYY-MM-DD hh:mm:ss"
     * @param t0            Initial time (sec)
     * @param t1            Final time (sec)
     * @param dt            Integration step (sec)
     */
    void integrate(Spacecraft& sc, const EnvironmentModel& env_model, const std::string& epoch, const Dimension::Time t0, const Dimension::Time t1, const Dimension::Time dt);
};


#endif //ORBITAL_PROPAGATOR_SIM_H
