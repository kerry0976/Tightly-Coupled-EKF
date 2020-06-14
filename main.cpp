
#include <stdint.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "GNSS_measurement.h"
#include "GNSS.h"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
using namespace Eigen;
using namespace std;

typedef Matrix<double, Dynamic, 7> MatrixX7d;
typedef Matrix<double, Dynamic, Dynamic> MatrixXd;

Matrix3d Skew(Vector3d w);
Matrix3f Skew(Vector3f w);
VectorXd EphemerisData2PosVelClock(GNSS_raw_measurement gnss_raw_measurement);
VectorXd EphemerisData2Satecef(float t,
                               uint32_t TOW, uint8_t L2, uint16_t week_No, uint8_t L2_Flag, uint8_t SV_Acc, uint8_t SV_Hlth,
                               double T_GD, uint16_t IODC, double t_OC, int8_t a_f2, double a_f1, double a_f0,
                               uint8_t IODE, double C_rs, double delta_n, double M_0, double C_uc, double ecc, double C_us,
                               double sqrt_A, double t_OE, double C_ic, double Omega_0, double C_is, double i_0, double C_rc,
                               double omega, double Omega_dot, double IDOT);

VectorXd GNSS_LS_pos_vel(MatrixXd gnss_measurement, int no_sat, Vector3d pEst_E_m_, Vector3d vEst_E_mps_);

int main(int arg_count, char *args[])
{

    GNSS_raw_measurement gnss_raw_measurement; 
    gnss_raw_measurement.AODO = 27900;         
    gnss_raw_measurement.Cic = -1.899898e-07;         
    gnss_raw_measurement.Cis = -1.192093e-07;            
    gnss_raw_measurement.Crc = 258.8438;         
    gnss_raw_measurement.Crs = 4.659375e+01;           
    gnss_raw_measurement.Cuc = 2.536923e-06;           
    gnss_raw_measurement.Cus = 4.950911e-06;         
    gnss_raw_measurement.FIT  = 0;             
    gnss_raw_measurement.IDOT = 7.35554e-11;          
    gnss_raw_measurement.IODC = 70;              
    gnss_raw_measurement.IODE = 70;               
    gnss_raw_measurement.L2  = 1;                   
    gnss_raw_measurement.L2P = 0;             
    gnss_raw_measurement.M0 = -8.16680097638e-01;            
    gnss_raw_measurement.Omega0 = 8.26951844152e-01;     
    gnss_raw_measurement.Omegad = -2.922775e-09;        
    gnss_raw_measurement.TOW17 =  81231;       
    gnss_raw_measurement.Tgd = -1.25729e-08;        
    gnss_raw_measurement.WN = 59;         
    gnss_raw_measurement.af0 = -2.6363600e-04;        
    gnss_raw_measurement.af1 = 1.068656e-11;      
    gnss_raw_measurement.af2 = 0; 
    gnss_raw_measurement.constellation = "GPS"; 
    gnss_raw_measurement.deltan = 1.876174e-09;        
    gnss_raw_measurement.doppler = -1252.290283;     
    gnss_raw_measurement.e = 0.015969;              
    gnss_raw_measurement.frame1 = "True";        
    gnss_raw_measurement.frame2  = "True";                
    gnss_raw_measurement.frame3 = "True";             
    gnss_raw_measurement.gnssid = 0;          
    gnss_raw_measurement.hlth = 0;              
    gnss_raw_measurement.i0 = 2.90582119487e-01;           
    gnss_raw_measurement.omega = 6.32357007824e-01;        
    gnss_raw_measurement.pseudorange = 22484920.970332;   
    gnss_raw_measurement.sqrtA = 5153.7164993;         
    gnss_raw_measurement.timestamp = 21600;    
    gnss_raw_measurement.toc  = 489600;          
    gnss_raw_measurement.toe  = 489600;          
    gnss_raw_measurement.ura = 0;         
    
    //GNSS gnss_measurement(gnss_raw_measurement);
    //VectorXd gnssProcessedmeasurement= gnss_measurement.getGNSS_measurement();
    VectorXd gnssProcessedmeasurement= EphemerisData2PosVelClock(gnss_raw_measurement);
    cout << endl;
    cout << "Pseudorange:  " << gnssProcessedmeasurement(0) << "m\n";
    cout << "Pseudorange rate: " << gnssProcessedmeasurement(1) << "m/s\n";
    cout << "x: " << gnssProcessedmeasurement(2) << "m\n";
    cout << "y: " << gnssProcessedmeasurement(3) << "m\n";
    cout << "z: " << gnssProcessedmeasurement(4) << "m\n";
    cout << "vx: " << gnssProcessedmeasurement(5) << "m\n";
    cout << "vy: " << gnssProcessedmeasurement(6) << "m\n";
    cout << "vz: " << gnssProcessedmeasurement(7) << "m\n";
   
    // Test Least Square Solution 
    vector<vector<string>> data;
    ifstream infile;
    infile.open("gnss_measurement.txt");
    while (infile)
    {
        string s;
        if (!getline(infile, s))
            break;

        istringstream ss(s);
        vector<string> record;

        while (ss)
        {
            string s;
            if (!getline(ss, s, ','))
                break;
            record.push_back(s);
        }

        data.push_back(record);
    }
    if (!infile.eof())
    {
        cerr << "Fooey!\n";
    }
    int no_row = data.size() - 1;
    int no_col = data[0].size();

    MatrixXd GNSS_measurement(data.size() - 1, 8);
    GNSS_measurement.setZero();

    for (int i = 0; i < no_row ; i++)
    {
        for (int j = 0; j < no_col ; j++)
        {
            // cout << stof(data[i][j]) * 1e+7 << " ";
            GNSS_measurement(i, j) = stod(data[i][j]) * 1e+7;
        }
    }
    cout << "\n";

    MatrixXd output(8,1);
    Vector3d pEst_E_m_, vEst_E_mps_;
    
    pEst_E_m_.setZero();
    vEst_E_mps_.setZero();
    cout << "Location of Vehicle using Linear Least Square: \n";
    cout << endl;
    output = GNSS_LS_pos_vel(GNSS_measurement, no_row, pEst_E_m_, vEst_E_mps_);
    cout << "x: " << output(0) << "(m), y: " <<output(1) << "(m), z: " <<  output(2) << "(m) " << endl; 
    cout << "vx: " <<  output(3) << "(m/s), vy: " <<output(4) << "(m/s), vz: " <<  output(5) << "(m/s) " << endl; 
    cout  << "clock offset: " << output(6) << "(m), clock shift: " <<output(7) << "(m/s) " << endl; 

    return 0;
}

VectorXd EphemerisData2Satecef(float t,
                               uint32_t TOW, uint8_t L2, uint16_t week_No, uint8_t L2_Flag, uint8_t SV_Acc, uint8_t SV_Hlth,
                               double T_GD, uint16_t IODC, double t_OC, int8_t a_f2, double a_f1, double a_f0,
                               uint8_t IODE, double C_rs, double delta_n, double M_0, double C_uc, double ecc, double C_us,
                               double sqrt_A, double t_OE, double C_ic, double Omega_0, double C_is, double i_0, double C_rc,
                               double omega, double Omega_dot, double IDOT)

{
    // All the equations are based ON: https://www.gps.gov/technical/icwg/IS-GPS-200H.pdf
    // pg. 104-105, Also Grove  p335-338

    // Process subframe 1,2,3 information
    double A_semiMajorAxis;        // Semi-major axis
    double n_0_computedMeanMotion; // Computed mean motion
    double n_correctedMeanMotion;  // Corrected mean motion
    double e_eccentricity;         // Eccentricity
    //double phi_k_argumentOfLattitude;   // Argument of latitude
    double M_0_trueAnomalyAtRef;
    double omega0_longitudeofAscendingNodeofOrbitPlane;
    double omega_argumentOfPerigee;
    double omegaDot_argumentOfPerigee;
    double i_0_inclinationAtRef;
    double iDot_rateOfInclination;

    A_semiMajorAxis = pow(sqrt_A, 2);
    n_0_computedMeanMotion = sqrt(MU / pow(A_semiMajorAxis, 3));
    n_correctedMeanMotion = n_0_computedMeanMotion + delta_n;
    e_eccentricity = ecc;
    M_0_trueAnomalyAtRef = M_0;
    omega0_longitudeofAscendingNodeofOrbitPlane = Omega_0;
    omega_argumentOfPerigee = omega;
    omegaDot_argumentOfPerigee = Omega_dot;
    i_0_inclinationAtRef = i_0;
    iDot_rateOfInclination = IDOT;

    // Compute the time from the ephemeris reference epoch
    double t_k_timeFromReferenceEpoch = t - t_OE;
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

    return pos_vel_Sat_ecef;
}

// Skew symmetric matrix from a given vector w
Matrix3d Skew(Vector3d w)
{
    Matrix3d C;

    C(0, 0) = 0.0;
    C(0, 1) = -w(2);
    C(0, 2) = w(1);
    C(1, 0) = w(2);
    C(1, 1) = 0.0;
    C(1, 2) = -w(0);
    C(2, 0) = -w(1);
    C(2, 1) = w(0);
    C(2, 2) = 0.0;

    return C;
}
Matrix3f Skew(Vector3f w)
{
    Matrix3f C;

    C(0, 0) = 0.0;
    C(0, 1) = -w(2);
    C(0, 2) = w(1);
    C(1, 0) = w(2);
    C(1, 1) = 0.0;
    C(1, 2) = -w(0);
    C(2, 0) = -w(1);
    C(2, 1) = w(0);
    C(2, 2) = 0.0;

    return C;
}

VectorXd GNSS_LS_pos_vel(MatrixXd gnss_measurement, int no_sat, Vector3d pEst_E_m_, Vector3d vEst_E_mps_)
{
    /*
   GNSS_LS_position_velocity - Calculates position, velocity, clock offset, 
   and clock drift using unweighted iterated least squares. Separate
   calculations are implemented for position and clock offset and for
   velocity and clock drift

    % Inputs:
    %   GNSS_measurements     GNSS measurement data:
    %     Column 0              Pseudo-range measurements (m)
    %     Column 1              Pseudo-range rate measurements (m/s)
    %     Columns 2-4           Satellite ECEF position (m)
    %     Columns 5-7           Satellite ECEF velocity (m/s)
    %   no_GNSS_meas          Number of satellites for which measurements are
    %                         supplied
    %   pEst_E_m_        prior predicted ECEF user position (m)
    %   vEst_E_mps_      prior predicted ECEF user velocity (m/s)
    %
    % Outputs:

        output is a 8 by 1 matrix containing the following: 

    %   est_r_ea_e            estimated ECEF user position (m)
    %   est_v_ea_e            estimated ECEF user velocity (m/s)
    %   est_clock             estimated receiver clock offset (m) and drift (m/s)
   */


    // Position and Clock OFFSET
    VectorXd x_pred(4, 1);
    x_pred.setZero();
    x_pred.segment(0, 3) = pEst_E_m_;
    x_pred(3) = vEst_E_mps_(0);

    // allocate space for estimation
    VectorXd delta_r(3, 1);
    delta_r.setZero();
    double approx_range = 0;
    double range = 0;
    double range_rate = 0;
    VectorXd pred_meas(no_sat, 1);
    pred_meas.setZero();
    MatrixXd H(no_sat, 4);
    H.setZero();
    Matrix<double, 3, 3> T_E2I;
    T_E2I.Identity();
    VectorXd u_as_E(3, 1);
    u_as_E.setZero();

    // output
    Matrix<double, 4, 1> x_est_1;
    x_est_1.setZero();
    Matrix<double, 4, 1> x_est_2;
    x_est_2.setZero();
    VectorXd output(8, 1);
    output.setZero();

    // set the flag for the while-loop
    double test_convergence = 1;
    while (test_convergence > 0.0001)
    {
        
        for (int j = 0; j < no_sat; j++)
        {
            Vector3d x_temp;
            x_temp(0) = gnss_measurement(j, 2);
            x_temp(1) = gnss_measurement(j, 3);
            x_temp(2) = gnss_measurement(j, 4);

            delta_r = x_temp - x_pred.segment(0, 3);
            approx_range = delta_r.norm();

            // Calculate frame rotation during signal transit time using (Grove2nd:8.36)
            T_E2I(0, 0) = 1;
            T_E2I(0, 1) = OMEGA_DOT_EARTH * approx_range / c;
            T_E2I(0, 2) = 0;
            T_E2I(1, 0) = -OMEGA_DOT_EARTH * approx_range / c;
            T_E2I(1, 1) = 1;
            T_E2I(1, 2) = 0;
            T_E2I(2, 0) = 0;
            T_E2I(2, 1) = 0;
            T_E2I(2, 2) = 1;

            delta_r = T_E2I.cast<double>() * x_temp - x_pred.segment(0, 3);
            range = delta_r.norm();
            pred_meas(j, 0) = range + x_pred(3);
            u_as_E = delta_r / range;
            H(j, 0) = -u_as_E(0);
            H(j, 1) = -u_as_E(1);
            H(j, 2) = -u_as_E(2);
            H(j, 3) = 1;
        }
        // least sqaure method
        x_est_1 = x_pred + (H.transpose() * H).inverse() * H.transpose() * (gnss_measurement.col(0) - pred_meas);
        test_convergence = (x_est_1 - x_pred).norm();
        x_pred = x_est_1;
    }
    
    // Earth rotation vector and matrix
    Vector3d omega_ie; // Earth rate vector
    Matrix3d Omega_ie; // Earth rate rotation matrix
    omega_ie(0) = 0.0;
    omega_ie(1) = 0.0;
    omega_ie(2) = OMEGA_DOT_EARTH;
    Omega_ie = Skew(omega_ie);

    // save the estimated postion and clock offset
    output(0) = x_est_1(0);
    output(1) = x_est_1(1);
    output(2) = x_est_1(2);
    output(6) = x_est_1(3); // record clock offset

    x_pred.setZero();
    pred_meas.setZero();
    x_pred.segment(0, 3) = vEst_E_mps_;
    x_pred(3) = 0;
    test_convergence = 1;
    u_as_E.setZero();
    
    while (test_convergence > 0.0001)
    {
        for (int j = 0; j < no_sat; j++)
        {
            Vector3d p_temp, v_temp;
            p_temp(0) = gnss_measurement(j, 2);
            p_temp(1) = gnss_measurement(j, 3);
            p_temp(2) = gnss_measurement(j, 4);
            v_temp(0) = gnss_measurement(j, 5);
            v_temp(1) = gnss_measurement(j, 6);
            v_temp(2) = gnss_measurement(j, 7);

            delta_r = p_temp - output.segment(0, 3);
            approx_range = delta_r.norm();

            // Calculate frame rotation during signal transit time using (Grove2nd:8.36)
            T_E2I(0, 0) = 1;
            T_E2I(0, 1) = OMEGA_DOT_EARTH * approx_range / c;
            T_E2I(0, 2) = 0;
            T_E2I(1, 0) = -OMEGA_DOT_EARTH * approx_range / c;
            T_E2I(1, 1) = 1;
            T_E2I(1, 2) = 0;
            T_E2I(2, 0) = 0;
            T_E2I(2, 1) = 0;
            T_E2I(2, 2) = 1;

            delta_r = T_E2I.cast<double>() * p_temp - output.segment(0, 3);
            range = delta_r.norm();
            u_as_E = delta_r / range;
     
            // Predict pseudo-range rate using (9.165)
            range_rate = u_as_E.transpose() * (T_E2I.cast<double>() * (v_temp + Omega_ie.cast<double>() * p_temp) -
                                               (x_pred.segment(0, 3) + Omega_ie.cast<double>() * output.segment(0, 3)));
       
            pred_meas(j, 0) = range_rate + x_pred(3);  
    
            H(j, 0) = -u_as_E(0);
            H(j, 1) = -u_as_E(1);
            H(j, 2) = -u_as_E(2);
            H(j, 3) = 1;
        }

        x_est_2 = x_pred + (H.transpose() * H).inverse() * H.transpose() * (gnss_measurement.col(1) - pred_meas);
        test_convergence = (x_est_2 - x_pred).norm();
        x_pred = x_est_2;
    }
    
    // save the estimated postion and clock offset
    output(3) = x_est_2(0);
    output(4) = x_est_2(1);
    output(5) = x_est_2(2);
    output(7) = x_est_2(3);
    return output;
}

// process raw measurement to give a 8 x 1 matrix (range, range rate, x,y,z,vx,vy,vz)
VectorXd EphemerisData2PosVelClock(GNSS_raw_measurement gnss_raw_measurement)
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

    VectorXd pos_vel_ecef_clock(8);
    double lambda = 2*c / (1575.4282e6);  // L1 according ublox8
    // https://www.u-blox.com/sites/default/files/products/documents/u-blox8-M8_ReceiverDescrProtSpec_%28UBX-13003221%29.pdf
    double PseudorangeRate = lambda * gnss_raw_measurement.doppler;
    pos_vel_ecef_clock << gnss_raw_measurement.pseudorange, PseudorangeRate, x, y, z, vx, vy, vz;
    // cout << x << y << y << endl;

    return pos_vel_ecef_clock;
}