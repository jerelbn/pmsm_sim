#pragma once

#include "common_cpp/common.h"
#include "common_cpp/logger.h"

class Motor {

enum StateIndices {
    ID, // direct current
    IQ, // quadrature current
    THETA, // angle
    OMEGA, // angular rate
    NS
};

enum InputIndices {
    VD, // direct voltage
    VQ, // quadrature voltage
    NI
};

public:
    Motor();
    Motor(const std::string &filename);
    ~Motor();

    void load(const std::string &filename);
    void update(const double &t);
    void computeControl(const double& t, const double& theta_c);

    const double* getInput() const { return input; }

private:

    void f(const double x[NS], const double u[NI], double dx[NS]);
    void log(const double &t);

    double t_prev; // previous timestamp

    // Electrical constants
    double L; // inductance
    double R; // resistance
    double Jm; // motor inertia
    double lam; // flux constant
    double np; // number of poles

    // Mechanical constants
    double Tsf; // static friction torque
    double mu; // viscous friction constant
    double Tl; // load torque
    double Jl; // load inertia

    // State array and change
    double state[NS];
    double dstate[NS];

    // PI controller
    double ctrl_update_rate;
    double input[NI], idq_c[NI];
    common::PID<double> pid_d;
    common::PID<double> pid_q;
    common::PID<double> pid_theta;

    common::Logger logger;
};
