//
// Created by Sylvain  on 3/31/20.
//

#include "Sim.h"

Sim::Sim(const std::string meta_kernel): m_meta_kernel(meta_kernel) {
    // Load SPICE kernels
    furnsh_c(m_meta_kernel.c_str());
}

Sim::~Sim() {
    // Clean up kernel space
    unload_c(m_meta_kernel.c_str());
}

void Sim::integrate(Spacecraft& sc, const EnvironmentModel& env_model, const std::string& epoch, const Dimension::Time t0, const Dimension::Time t1, const Dimension::Time dt) {
    SpiceDouble et;
    str2et_c(epoch.c_str(), &et);
    EqMotion eom(sc, env_model, et);
    Integrator::integrate(sc, eom, t0.data(), t1.data(), dt.data());
}