//
// Created by Sylvain  on 4/7/20.
//

#ifndef ORBITAL_PROPAGATOR_SENSOR_H
#define ORBITAL_PROPAGATOR_SENSOR_H

#include <Eigen/Dense>
#include "State.h"
#include "EnvironmentModel.h"

class EnvironmentModel;

/** Base class representing a sensor */
class Sensor {
public:
    Sensor() {}
    virtual ~Sensor() {}

    /**
     * @param state         State vector of the spacecraft
     * @param et            Ephemeris time
     * @param env_model     Reference to the environment model
     * @return              Measurement in the sensor-fixed frame.
     */
    virtual Eigen::Vector3d measure(const State& state, double et, const EnvironmentModel& env_model) = 0;
};


#endif //ORBITAL_PROPAGATOR_SENSOR_H
