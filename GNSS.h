#ifndef GNSSm
#define GNSSm

#include <stdint.h>
#include <math.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/LU>

using namespace Eigen;

#include "GNSS_measurement.h"


const double  MU = 3.986005e14; //m^3 / sec^2
const double OMEGA_DOT_EARTH = 7.2921151467e-5; // rad/sec
const double c = 299792458; // Speed of light in m/s

class GNSS {
public:
    GNSS()  {};
    GNSS(const struct GNSS_raw_measurement&);
    VectorXd getGNSS_measurement();
    VectorXd EphemerisData2Satecef(GNSS_raw_measurement gnss_raw_measurement);
    // ~GNSS() {};

private:
    Matrix<double,8,1> gnss_measurement; // processed measurement 8 x 1 matrix. Each row: pesudorange, range-rate, position, velocity
    VectorXd pos_vel_Sat_ecef; // position and velocity in ECEF
    double Pseudorange;
    double PseudorangeRate;
};

#endif