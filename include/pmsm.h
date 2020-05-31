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
    void update(const double &t, const double u[NI]);

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

    // State array and derivative
    double state[NS];
    double state_dot[NS];

    common::Logger logger;
};
