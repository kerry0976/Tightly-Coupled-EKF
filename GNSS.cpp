#include "GNSS_measurement.h"
#include "GNSS.h"

#include <stdint.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/LU>

using namespace std;
using namespace Eigen;

GNSS::GNSS(const GNSS_raw_measurement &gnss_raw_measurement)
{
    // cout << "contructor\n";
    pos_vel_Sat_ecef = EphemerisData2Satecef(gnss_raw_measurement);
    Pseudorange = gnss_raw_measurement.pseudorange;
    double lambda = 2*c / (1575.4282e6);  // L1 according ublox8
    // https://www.u-blox.com/sites/default/files/products/documents/u-blox8-M8_ReceiverDescrProtSpec_%28UBX-13003221%29.pdf
    PseudorangeRate = lambda * gnss_raw_measurement.doppler;
}

VectorXd GNSS::getGNSS_measurement()
{
   
    gnss_measurement(0) = Pseudorange;
    gnss_measurement(1) = PseudorangeRate;
    gnss_measurement(2) = pos_vel_Sat_ecef(0);
    gnss_measurement(3) = pos_vel_Sat_ecef(1);
    gnss_measurement(4) = pos_vel_Sat_ecef(2);
    gnss_measurement(5) = pos_vel_Sat_ecef(3);
    gnss_measurement(6) = pos_vel_Sat_ecef(4);
    gnss_measurement(7) = pos_vel_Sat_ecef(5);
    return gnss_measurement;
}

// Ephemeris Data (subframe1,2,3) to Satellite ecef x, y, z in meter, vx, vy, vz in m/s
VectorXd GNSS::EphemerisData2Satecef(GNSS_raw_measurement gnss_raw_measurement)
{
    // All the equations are based ON: https://www.gps.gov/technical/icwg/IS-GPS-200H.pdf
    // pg. 104-105, Also Grove  p335-338

    // Process subframe 1,2,3 information
    double A_semiMajorAxis = pow(gnss_raw_measurement.sqrtA, 2);        // Semi-major axis
    double n_0_computedMeanMotion = sqrt(MU / pow(A_semiMajorAxis, 3)); // Computed mean motion
    double n_correctedMeanMotion = n_0_computedMeanMotion + gnss_raw_measurement.deltan;// Corrected mean motion
    double e_eccentricity = gnss_raw_measurement.e; // Eccentricity
    //double phi_k_argumentOfLattitude;   // Argument of latitude
    double M_0_trueAnomalyAtRef = gnss_raw_measurement.M0;
    double omega0_longitudeofAscendingNodeofOrbitPlane = gnss_raw_measurement.Omega0;
    double omega_argumentOfPerigee = gnss_raw_measurement.omega;
    double omegaDot_argumentOfPerigee = gnss_raw_measurement.Omegad;
    double i_0_inclinationAtRef = gnss_raw_measurement.i0;
    double iDot_rateOfInclination = gnss_raw_measurement.IDOT;

    double t_OE = gnss_raw_measurement.toe;
    // 2nd harmonic terms
    double C_us = gnss_raw_measurement.Cus;
    double C_uc = gnss_raw_measurement.Cuc;
    double C_rs = gnss_raw_measurement.Crs;
    double C_rc = gnss_raw_measurement.Crc;
    double C_is = gnss_raw_measurement.Cis;
    double C_ic = gnss_raw_measurement.Cic;

    // Compute the time from the ephemeris reference epoch
    double t_k_timeFromReferenceEpoch = gnss_raw_measurement.timestamp - t_OE;
    // Correct that time for end-of-week crossovers
    if (t_k_timeFromReferenceEpoch > 302400)
    {
        t_k_timeFromReferenceEpoch -= 604800;
    }
    if (t_k_timeFromReferenceEpoch < -302400)
    {
        t_k_timeFromReferenceEpoch += 604800;
    }

    // Compute the mean anomaly
    double M_k_meanAnomaly = M_0_trueAnomalyAtRef + n_correctedMeanMotion * t_k_timeFromReferenceEpoch;

    // Below, we iteratively solve for E_k_eccentricAnomaly using Newton-Raphson method
    double solutionError = 1000000.;
    double E_k_eccentricAnomaly = 1.;
    double currentDerivative = 0;
    int iterationCount = 0;

    solutionError = (E_k_eccentricAnomaly -
                     (e_eccentricity * sin(E_k_eccentricAnomaly)) -
                     M_k_meanAnomaly);

    while ((fabs(solutionError) > 1.0e-6) &&
           iterationCount < 1000)
    {
        currentDerivative = (1.0 - (e_eccentricity * cos(E_k_eccentricAnomaly)));
        E_k_eccentricAnomaly = E_k_eccentricAnomaly - solutionError / currentDerivative;

        solutionError = (E_k_eccentricAnomaly -
                         (e_eccentricity * sin(E_k_eccentricAnomaly)) -
                         M_k_meanAnomaly);
        iterationCount += 1;
        //   if (VERBOSE)
        //   {
        //     std::cout<< "Iteration #: " << iterationCount << " Error: " << solutionError << std::endl;
        //   }
    }
    double cos_E_k = cos(E_k_eccentricAnomaly);
    double sin_E_k = sin(E_k_eccentricAnomaly);
    double nu_k_trueAnomaly = atan2(
        (sqrt(1.0 - pow(e_eccentricity, 2)) * sin_E_k) /
            (1.0 - (e_eccentricity * cos_E_k)),
        (cos_E_k - e_eccentricity) /
            (1.0 - e_eccentricity * cos_E_k));

    double phi_k_argumentOfLatitude = nu_k_trueAnomaly + omega_argumentOfPerigee;

    // Compute the corrective 2nd order terms
    double sin2PhiK = sin(2.0 * phi_k_argumentOfLatitude);
    double cos2PhiK = cos(2.0 * phi_k_argumentOfLatitude);

    double deltaU_argumentOfLatCorrection = (C_us * sin2PhiK) + (C_uc * cos2PhiK);
    double deltaR_radiusCorrection = (C_rs * sin2PhiK) + (C_rc * cos2PhiK);
    double deltaI_inclinationCorrection = (C_is * sin2PhiK) + (C_ic * cos2PhiK);

    // Now compute the updated corrected orbital elements
    double u_argumentOfLat = phi_k_argumentOfLatitude + deltaU_argumentOfLatCorrection;
    double r_radius = (A_semiMajorAxis * (1.0 - (e_eccentricity * cos_E_k))) + deltaR_radiusCorrection;
    double i_inclination =
        i_0_inclinationAtRef +
        (iDot_rateOfInclination * t_k_timeFromReferenceEpoch) +
        deltaI_inclinationCorrection;

    // Compute the satellite position within the orbital plane
    double xPositionOrbitalPlane = r_radius * cos(u_argumentOfLat);
    double yPositionOrbitalPlane = r_radius * sin(u_argumentOfLat);
    double omegaK_longitudeAscendingNode =
        omega0_longitudeofAscendingNodeofOrbitPlane +
        ((omegaDot_argumentOfPerigee - OMEGA_DOT_EARTH) * t_k_timeFromReferenceEpoch) -
        (OMEGA_DOT_EARTH * t_OE);

    double sinOmegaK = sin(omegaK_longitudeAscendingNode);
    double cosOmegaK = cos(omegaK_longitudeAscendingNode);

    double sinIK = sin(i_inclination);
    double cosIK = cos(i_inclination);
    // Earth-fixed coordinates:
    double x = (xPositionOrbitalPlane * cosOmegaK) - (yPositionOrbitalPlane * cosIK * sinOmegaK);
    double y = (xPositionOrbitalPlane * sinOmegaK) + (yPositionOrbitalPlane * cosIK * cosOmegaK);
    double z = (yPositionOrbitalPlane * sinIK);

    // ECEF velocity calculation:
    double E_dot_k_eccentricAnomaly = n_correctedMeanMotion / (1.0 - (e_eccentricity * cos_E_k));     // Eq.(8.21)
    double phi_dot_k_argumentOfLatitude = sin(nu_k_trueAnomaly) / sin_E_k * E_dot_k_eccentricAnomaly; // Eq.(8.22)
    double r_dot_o_os = (A_semiMajorAxis * e_eccentricity * sin_E_k) * E_dot_k_eccentricAnomaly +
                        2 * ((C_rs * cos2PhiK) - (C_rc * sin2PhiK)) * phi_dot_k_argumentOfLatitude;     // Eq.(8.23a)
    double u_dot_o_os = (1 + 2 * C_us * cos2PhiK - 2 * C_uc * sin2PhiK) * phi_dot_k_argumentOfLatitude; // Eq.(8.23b)

    double x_dot_o_os = r_dot_o_os * cos(u_argumentOfLat) - r_radius * u_dot_o_os * sin(u_argumentOfLat); // Eq.(8.24a)
    double y_dot_o_os = r_dot_o_os * sin(u_argumentOfLat) + r_radius * u_dot_o_os * cos(u_argumentOfLat); // Eq.(8.24b)

    double omega_dot_K_longitudeAscendingNode = omegaDot_argumentOfPerigee - OMEGA_DOT_EARTH;                                       // Eq. (8.25)
    double i_dot_inclination = iDot_rateOfInclination + 2 * ((C_is * cos2PhiK) - (C_ic * sin2PhiK)) * phi_dot_k_argumentOfLatitude; // Eq. (8.26)

    // Eq. (8.27)
    double vx = (x_dot_o_os * cosOmegaK - y_dot_o_os * cosIK * sinOmegaK + i_dot_inclination * yPositionOrbitalPlane * sinIK * sinOmegaK) -
                omega_dot_K_longitudeAscendingNode * (xPositionOrbitalPlane * sinOmegaK + yPositionOrbitalPlane * cosIK * cosOmegaK);
    double vy = (x_dot_o_os * sinOmegaK + y_dot_o_os * cosIK * cosOmegaK - i_dot_inclination * yPositionOrbitalPlane * sinIK * cosOmegaK) -
                omega_dot_K_longitudeAscendingNode * (-xPositionOrbitalPlane * cosOmegaK + yPositionOrbitalPlane * cosIK * sinOmegaK);
    double vz = (y_dot_o_os * sinIK + i_dot_inclination * yPositionOrbitalPlane * cosIK);

    VectorXd pos_vel_Sat_ecef(6);
    pos_vel_Sat_ecef << x, y, z, vx, vy, vz;
    // cout << x << y << y << endl;
    cout << "Emphermris function end\n";
    return pos_vel_Sat_ecef;
}