#include "common_cpp/common.h"
#include "common_cpp/progress_bar.h"
#include "pmsm.h"

int main()
{
    int seed;
    double t(0), tf, dt;
    common::getYamlNode("seed", std::string(PARAM_DIR)+"/pmsm.yaml", seed);
    common::getYamlNode("tf", std::string(PARAM_DIR)+"/pmsm.yaml", tf);
    common::getYamlNode("dt", std::string(PARAM_DIR)+"/pmsm.yaml", dt);
    if (seed < 0)
        seed = time(NULL);
    std::srand(seed);

    // Create progress bar
    common::ProgressBar prog_bar;
    prog_bar.init(tf / dt, 40);

    // Create vehicles, controllers, estimators, sensor packages
    Motor motor(std::string(PARAM_DIR)+"/pmsm.yaml");

    // Main simulation loop
    while (t <= tf+dt)
    {
        // Update motor control, state, etc.
        double u[] = {0.0, 1.0};
        motor.update(t, u);

        // Update time step
        t += dt;
        prog_bar.print(t / dt);
    }
    prog_bar.finished();

    return 0;
}