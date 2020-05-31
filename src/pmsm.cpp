#include "pmsm.h"

Motor::Motor() : t_prev{},L{},R{},Jm{},lam{},np{},Tsf{},mu{},Tl{},Jl{},state{} {}

Motor::Motor(const std::string &filename) : t_prev{},L{},R{},Jm{},lam{},np{},Tsf{},mu{},Tl{},Jl{},state{} {
    load(filename);
}

Motor::~Motor() {}

void Motor::load(const std::string &filename)
{
    // Load all parameters
    common::getYamlNode("mc", filename, mc_);
    common::getYamlNode("mp", filename, mp_);
    common::getYamlNode("Jy", filename, Jy_);
    common::getYamlNode("Jz", filename, Jz_);
    common::getYamlNode("Jp", filename, Jp_);
    common::getYamlNode("Jm", filename, Jm_);
    common::getYamlNode("Km", filename, Km_);
    common::getYamlNode("Lm", filename, Lm_);
    common::getYamlNode("Rm", filename, Rm_);
    common::getYamlNode("L", filename, L_);
    common::getYamlNode("l", filename, l_);
    common::getYamlNode("r", filename, r_);
    common::getYamlNode("bm", filename, bm_);
    common::getYamlNode("bp", filename, bp_);

    // Initialize loggers and log initial data
    std::string logname_true_state;
    common::getYamlNode("logname_true_state", filename, logname_true_state);
    state_log_.open(logname_true_state);
}

void Motor::update(const double &t, const double u[NI])
{
    double dt = t - t_prev;
    t_prev = t;

    if (t > 0)
    {
        // 4th order Runge-Kutta integration
        common::rk4<double, NUM_STATES, NUM_INPUTS>(std::bind(&Motor::f, this,
                                                              std::placeholders::_1,
                                                              std::placeholders::_2,
                                                              std::placeholders::_3),
                                                    dt, x_, u, dx_);
        x_ += dx_;
    }

    log(t);
}

void Motor::f(const xVector &x, const uVector &u, xVector &dx)
{
    // Constants
    static double g = common::gravity;
    static double c0 = mp_*l_;
    static double c1 = Jy_ + c0*l_;
    static double c2 = c0*r_;
    static double c22 = c2*c2;
    static double c3 = 2.0*Jm_*c1;
    static double c4 = 2.0*Jm_*c2;
    static double c5 = c0*g;

    // Unpack states/inputs for readability
    double theta = x_(THETA);
    double dtheta = x_(DTHETA);
    double omegal = x_(OMEGAL);
    double omegar = x_(OMEGAR);
    double ql = x_(QL);
    double qr = x_(QR);
    double Vl = u(VL);
    double Vr = u(VR);

    // Common calcs
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double cos2_theta = cos_theta*cos_theta;
    double dtheta2 = dtheta*dtheta;
    double v0 = c22*cos2_theta;
    double v1 = (c1*dtheta2 - c5)*sin_theta - bp_*dtheta;
    double v2 = c4*cos_theta*v1;
    double v3 = Km_*(ql + qr);
    double v4 = bm_*(omegal + omegar);
    double v5 = Jm_*c2*cos_theta;
    double v6 = c2*cos_theta;
    double denom = Jm_*(Jm_*c1 - v0);
    
    // Equations of motion
    dx(DX) = -(r_/2.0)*((c3 + v0*(Jm_ - 1.0))*(v3 - v4) + v2)/denom;
    dx(DPSI) = (r_/L_)*((c3 - v0*(Jm_ + 1.0))*(Km_*(qr - ql) + bm_*(omegal - omegar)))/denom;
    dx(THETA) = dtheta;
    dx(DTHETA) = (c5*sin_theta-bp_*dtheta)/c1 - v6/(2.0*c1)*((c3 + v0*(Jm_-1.0))*(v3 - v4) + v2)/denom;
    dx(OMEGAL) = ((c3 - v0)*(Km_*ql - bm_*omegal) + v5*(v6*(Km_*qr - bm_*omegar) + v1))/denom;
    dx(OMEGAR) = ((c3 - v0)*(Km_*qr - bm_*omegar) + v5*(v6*(Km_*ql - bm_*omegal) + v1))/denom;
    dx(QL) = (Vl - Rm_*ql - Km_*omegal)/Lm_;
    dx(QR) = (Vr - Rm_*qr - Km_*omegar)/Lm_;
}

void Motor::log(const double &t)
{
    state_log_.log(t);
    state_log_.logMatrix(x_);
}