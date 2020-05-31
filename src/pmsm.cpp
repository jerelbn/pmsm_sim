#include "pmsm.h"

Motor::Motor() : t_prev{},L{},R{},Jm{},lam{},np{},Tsf{},mu{},Tl{},Jl{},state{},dstate{},input{},idq_c{} {}

Motor::Motor(const std::string &filename) : t_prev{},L{},R{},Jm{},lam{},np{},Tsf{},mu{},Tl{},Jl{},state{},dstate{},input{},idq_c{} {
    load(filename);
}

Motor::~Motor() {}

void Motor::load(const std::string &filename)
{
    // Load all parameters
    common::getYamlNode("motor_inductance", filename, L);
    common::getYamlNode("motor_resistance", filename, R);
    common::getYamlNode("motor_inertia", filename, Jm);
    common::getYamlNode("motor_flux_constant", filename, lam);
    common::getYamlNode("motor_poles", filename, np);
    common::getYamlNode("motor_static_friction_torque", filename, Tsf);
    common::getYamlNode("motor_viscous_friction_constant", filename, mu);
    common::getYamlNode("motor_load_torque", filename, Tl);
    common::getYamlNode("motor_load_inertia", filename, Jl);

    common::getYamlNode("ctrl_update_rate", filename, ctrl_update_rate);

    pid_d.init(500.0, 10.0, 0.0, 12.0, -12.0, 0.2, -0.2, 0.5);
    pid_q.init(50.0, 0.0, 0.0, 12.0, -12.0, 0.2, -0.2, 0.5);
    pid_theta.init(1.5, 0.01, 0.0, 12.0, -12.0, 0.1, -0.1, 0.5);

    // Initialize loggers and log initial data
    std::string logname_true_state;
    common::getYamlNode("logname_true_state", filename, logname_true_state);
    logger.open(logname_true_state);
}

void Motor::update(const double &t)
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
                                        dt, state, input, dstate);
        state[ID] += dstate[ID];
        state[IQ] += dstate[IQ];
        state[THETA] += dstate[THETA];
        state[OMEGA] += dstate[OMEGA];
    }

    log(t);
}

void Motor::computeControl(const double& t, const double& theta_c)
{
    static double ctrl_t_prev = 0;
    double dt = common::roundDecimal(t - ctrl_t_prev, 8);

    if (t == 0 || dt >= 1.0 / ctrl_update_rate)
    {
        idq_c[1] = pid_theta.run(dt, state[THETA], theta_c, true);
        input[VD] = pid_d.run(dt, state[ID], idq_c[0], true);
        input[VQ] = pid_q.run(dt, state[IQ], idq_c[1], true);
        ctrl_t_prev = t;
    }
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
    logger.log(t, state[ID], state[IQ], state[THETA], state[OMEGA], input[VD], input[VQ], idq_c[0], idq_c[1]);
}