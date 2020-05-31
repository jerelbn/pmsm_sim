#include "pmsm.h"

Motor::Motor() : t_prev{},L{},R{},Jm{},lam{},np{},Tsf{},mu{},Tl{},Jl{},state{},state_dot{} {}

Motor::Motor(const std::string &filename) : t_prev{},L{},R{},Jm{},lam{},np{},Tsf{},mu{},Tl{},Jl{},state{} {
    load(filename);
}

Motor::~Motor() {}

void Motor::load(const std::string &filename)
{
    // Load all parameters
    common::getYamlNode("motor_inductance", filename, L);
    common::getYamlNode("motor_resitance", filename, R);
    common::getYamlNode("motor_inertia", filename, Jm);
    common::getYamlNode("motor_flux_constant", filename, lam);
    common::getYamlNode("motor_poles", filename, np);
    common::getYamlNode("motor_static_friction_torque", filename, Tsf);
    common::getYamlNode("motor_viscous_friction_constant", filename, mu);
    common::getYamlNode("motor_load_torque", filename, Tl);
    common::getYamlNode("motor_load_inertia", filename, Jl);

    // Initialize loggers and log initial data
    std::string logname_true_state;
    common::getYamlNode("logname_true_state", filename, logname_true_state);
    logger.open(logname_true_state);
}

void Motor::update(const double &t, const double u[NI])
{
    double dt = t - t_prev;
    t_prev = t;

    if (t > 0)
    {
        // 4th order Runge-Kutta integration
        common::rk4<double, NS, NI>(std::bind(&Motor::f, this,
                                        std::placeholders::_1,
                                        std::placeholders::_2,
                                        std::placeholders::_3),
                                      dt, state, u, state_dot);
        state[ID] += state_dot[ID];
        state[IQ] += state_dot[IQ];
        state[THETA] += state_dot[THETA];
        state[OMEGA] += state_dot[OMEGA];
    }

    log(t);
}

// https://www.mathworks.com/help/physmod/sps/powersys/ref/permanentmagnetsynchronousmachine.html
void Motor::f(const double x[NS], const double u[NI], double dx[NS])
{
    double Te = 1.5*np*lam*x[IQ]; // Electrical torque
    dx[ID] = (u[VD] - R*x[ID])/L + np*x[OMEGA]*x[IQ];
    dx[IQ] = (u[VQ] - R*x[IQ])/L - (lam/L + x[ID])*np*x[OMEGA];
    dx[THETA] = x[OMEGA];
    dx[OMEGA] = (Te - Tsf - mu*x[OMEGA] - Tl)/Jl;
}

void Motor::log(const double &t)
{
    logger.log(t, state[ID], state[IQ], state[THETA], state[OMEGA]);
}