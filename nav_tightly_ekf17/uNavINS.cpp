/*
Copyright (c) 2016 - 2020 Regents of the University of Minnesota and Bolder Flight Systems Inc.
MIT License; See LICENSE.md for complete details

Adapted for 17-state Tightly-coupled EKF: Kerry Sun

Adapted from prior versions
Copyright 2011 Regents of the University of Minnesota. All rights reserved.
Original Author: Adhika Lie, Gokhan Inalhan, Demoz Gebre, Jung Soon Jang

Reference Frames and Coordinates from nav-functions()
I - ECI (Earch Center Inertial): origin at Earth center
E - ECEF (Earch Center Earth Fixed): origin at Earth center
D - Geodetic: origin at Earth center, Uses earth ellisoid definition (example WGS84)
G - Geocentric: origin at Earth center, Uses spheroid definition
L - Local Level: origin at specified reference, [x- North, y- East, z- Down]
B - Body: origin at Body CG, [x- Fwd, y- Starboard, z- Down]

All units meters and radians
"Acceleration" is actually "specific gravity", ie. gravity is removed.
*/
#include <iostream>
#include "uNavINS.h"
using namespace std; 

void uNavINS::Configure()
{
  // Observation matrix (H)
  //H_.setZero();
  //H_.block(0, 0, 5, 5) = I5;

  // Covariance of the Process Noise (associated with TimeUpdate())
  Rw_.setZero();
  Rw_.block(0, 0, 3, 3) = (aNoiseSigma_mps2 * aNoiseSigma_mps2) * I3;
  Rw_.block(3, 3, 3, 3) = (wNoiseSigma_rps * wNoiseSigma_rps) * I3;
  Rw_.block(6, 6, 3, 3) = 2.0f * (aMarkovSigma_mps2 * aMarkovSigma_mps2) / aMarkovTau_s * I3;
  Rw_.block(9, 9, 3, 3) = 2.0f * (wMarkovSigma_rps * wMarkovSigma_rps) / wMarkovTau_s * I3;
  Rw_(12,12) = 1.0f * pseudorangeNoiseSigma_m * pseudorangeNoiseSigma_m; 
  Rw_(13,13) = 1.0f * pseudorangeRateNoiseSigma_mps * pseudorangeRateNoiseSigma_mps;

  // Covariance of the Observation Noise (associated with MeasUpdate())
  // R_.setZero();
  // R_.block(0, 0, 2, 2) = (pNoiseSigma_NE_m * pNoiseSigma_NE_m) * I2;
  // R_(2, 2) = (pNoiseSigma_D_m * pNoiseSigma_D_m);
  // R_.block(3, 3, 2, 2) = (vNoiseSigma_NE_mps * vNoiseSigma_NE_mps) * I2;
  // R_(5, 5) = (vNoiseSigma_D_mps * vNoiseSigma_D_mps);

  // Initial Innovation Covariance Estimate (S)
  // S_.setZero();

  // Initial Covariance Estimate (P)
  P_.setZero();
  P_.block(0, 0, 3, 3) = (pErrSigma_Init_m * pErrSigma_Init_m) * I3;
  P_.block(3, 3, 3, 3) = (vErrSigma_Init_mps * vErrSigma_Init_mps) * I3;
  P_.block(6, 6, 2, 2) = (attErrSigma_Init_rad * attErrSigma_Init_rad) * I2;
  P_(8, 8) = (hdgErrSigma_Init_rad * hdgErrSigma_Init_rad);
  P_.block(9, 9, 3, 3) = (aBiasSigma_Init_mps2 * aBiasSigma_Init_mps2) * I3;
  P_.block(12, 12, 3, 3) = (wBiasSigma_Init_rps * wBiasSigma_Init_rps) * I3;
  P_(15, 15) = clockSigma_Init_m * clockSigma_Init_m;
  P_(16, 16) = clockrateSigma_Init_mps *  clockrateSigma_Init_mps;
}

void uNavINS::Initialize(Vector3f wMeas_B_rps, Vector3f aMeas_B_mps2, 
                         Vector3f magMeas_B_uT,Vector3d pMeas_D_rrm,
                         Vector3f vMeas_L_mps)
{

  // Initialize Position and Velocity
  //MatrixXd output(8,1);
  //output = GNSS_LS_pos_vel(GNSS_measurement, no_row, pMeas_D_rrm, vMeas_L_mps);
  // Initialize GNSS clock offset/phase shift
  // This should be initialized properly with a linear least square, todo 
  clockBias_m_ = 1.0002e+04; // clock bias
  clockRateBias_mps_ = 99.9828; // clock rate bias
  pEst_D_rrm_ = pMeas_D_rrm; // Position in LLA (rad, rad, m)
  vEst_L_mps_ = vMeas_L_mps; // Velocity in NED

  // Initialize sensor biases
  wBias_rps_ = wMeas_B_rps;
  aBias_mps2_.setZero();

  // New Specific forces and Rotation Rate
  aEst_B_mps2_ = aMeas_B_mps2 - aBias_mps2_;
  wEst_B_rps_ = wMeas_B_rps - wBias_rps_;

  // Initial attitude, roll and pitch
  Vector3f aEst_B_nd = aEst_B_mps2_ / aEst_B_mps2_.norm(); // Normalize to remove the 1g sensitivity
  euler_BL_rad_(1) = asinf(aEst_B_nd(0));
  euler_BL_rad_(0) = -asinf(aEst_B_nd(1) / cosf(euler_BL_rad_(1)));
  



  // Magnetic horizontal Y (right side) and X (forward) corrections due to roll and pitch angle
  Vector3f magMeas_B_nd = magMeas_B_uT / magMeas_B_uT.norm();
  float magY = magMeas_B_nd(1) * cosf(euler_BL_rad_(0)) - magMeas_B_nd(2) * sinf(euler_BL_rad_(0));
  float magX = magMeas_B_nd(0) * cosf(euler_BL_rad_(1)) + (magMeas_B_nd(1) * sinf(euler_BL_rad_(0)) + magMeas_B_nd(2) * cosf(euler_BL_rad_(0))) * sinf(euler_BL_rad_(1));

  // Estimate initial heading
  euler_BL_rad_(2) = -atan2f(magY, magX);

  // Euler to quaternion
  quat_BL_ = Euler2Quat(euler_BL_rad_);

  // set initialized flag
  initialized_ = true;
}

//void uNavINS::Update(uint64_t t_us, unsigned long timeWeek, Vector3f wMeas_B_rps, Vector3f aMeas_B_mps2, Vector3f ////magMeas_B_uT, Vector3d pMeas_D_rrm, Vector3f vMeas_L_mps)
void uNavINS::Update(uint64_t t_us, unsigned long timeWeek, 
                     Vector3f wMeas_B_rps, Vector3f aMeas_B_mps2, 
                     Vector3f magMeas_B_uT, GNSS_raw_measurement gnss_raw_measurement)
{
  // change in time
  dt_s_ = ((float)(t_us - tPrev_us_)) / 1e6;
  tPrev_us_ = t_us;

  // Catch large dt
  if (dt_s_ > 0.1)
  {
    dt_s_ = 0.1;
  }

  // A-priori accel and rotation rate estimate
  aEst_B_mps2_ = aMeas_B_mps2 - aBias_mps2_;
  wEst_B_rps_ = wMeas_B_rps - wBias_rps_;

  // Kalman Time Update (Prediction)
  TimeUpdate();

  // Gps measurement update, if TOW increased
  if ((timeWeek - timeWeekPrev_) > 0)
  {
    timeWeekPrev_ = timeWeek;

    // Process raw GNSS measurements (from ephemeris to a no_meas by 8 matrix )
    MatrixXd gnss_measurement;
    gnss_measurement = getGNSSmeasurement(gnss_raw_measurement, timeWeek);

    // Kalman Measurement Update
    MeasUpdate17(gnss_measurement);

    // Kalman Measurement Update
    //MeasUpdate(pMeas_D_rrm, vMeas_L_mps);

    // Post-priori accel and rotation rate estimate, biases updated in MeasUpdate()
    aEst_B_mps2_ = aMeas_B_mps2 - aBias_mps2_;
    wEst_B_rps_ = wMeas_B_rps - wBias_rps_;
  
  }

  // Euler angles from quaternion
  euler_BL_rad_ = Quat2Euler(quat_BL_);
}

void uNavINS::TimeUpdate()
{
  // Attitude Update
  Quaternionf dQuat_BL = Quaternionf(1.0, 0.5f * wEst_B_rps_(0) * dt_s_, 0.5f * wEst_B_rps_(1) * dt_s_, 0.5f * wEst_B_rps_(2) * dt_s_);
  quat_BL_ = (quat_BL_ * dQuat_BL).normalized();

  // Avoid quaternion flips sign
  if (quat_BL_.w() < 0)
  {
    quat_BL_ = Quaternionf(-quat_BL_.w(), -quat_BL_.x(), -quat_BL_.y(), -quat_BL_.z());
  }

  // Compute DCM (Body to/from NED) Transformations from Quaternion
  Matrix3f T_L2B = Quat2DCM(quat_BL_);
  Matrix3f T_B2L = T_L2B.transpose();

  // Velocity Update
  Vector3f aGrav_mps2 = Vector3f(0.0, 0.0, G);
  vEst_L_mps_ += dt_s_ * (T_B2L * aEst_B_mps2_ + aGrav_mps2);

  // Position Update
  Vector3f pDot_D = L2D_Rate(vEst_L_mps_, pEst_D_rrm_);
  pEst_D_rrm_ += (dt_s_ * pDot_D).cast<double>();

  // Clock Offset Update
  clockBias_m_ += clockRateBias_mps_ * dt_s_; 

  // Assemble the Jacobian (state update matrix)
  Matrix<float, 17, 17> Fs;
  Fs.setZero();
  Fs.block(0, 3, 3, 3) = I3;
  Fs(5, 2) = -2.0f * G / EARTH_RADIUS;
  Fs.block(3, 6, 3, 3) = -2.0f * T_B2L * Skew(aEst_B_mps2_);
  Fs.block(3, 9, 3, 3) = -T_B2L;
  Fs.block(6, 6, 3, 3) = -Skew(wEst_B_rps_);
  Fs.block(6, 12, 3, 3) = -0.5f * I3;
  Fs.block(9, 9, 3, 3) = -1.0f / aMarkovTau_s * I3;   // ... Accel Markov Bias
  Fs.block(12, 12, 3, 3) = -1.0f / wMarkovTau_s * I3; // ... Rotation Rate Markov Bias
  Fs(16,15) = 1.0f;

  // State Transition Matrix
  Matrix<float, 17, 17> PHI = I17 + Fs * dt_s_;

  // Process Noise Covariance (Discrete approximation)
  Matrix<float, 17, 14> Gs;
  Gs.setZero();
  Gs.block(3, 0, 3, 3) = -T_B2L;
  Gs.block(6, 3, 3, 3) = -0.5f * I3;
  Gs.block(9, 6, 3, 3) = I3;
  Gs.block(12, 9, 3, 3) = I3;
  Gs.block(15, 12, 2, 2) = I2;

  // Discrete Process Noise
  Matrix<float, 17, 17> Q;
  Q.setZero();
  Q = PHI * dt_s_ * Gs * Rw_ * Gs.transpose();
  Q = 0.5f * (Q + Q.transpose());

  // Covariance Time Update
  P_ = PHI * P_ * PHI.transpose() + Q;
  P_ = 0.5f * (P_ + P_.transpose());
}

// process raw gnss measurement
MatrixXd uNavINS::getGNSSmeasurement(GNSS_raw_measurement gnss_raw_measurement, unsigned long timeWeek)
{
  // pre-allocate matrix for the processed gnss measurement
  //int no_sat= gnss_raw_measurement.rows();
  int no_sat= 10;
  MatrixXd gnss_measurement(no_sat, 8);
  VectorXd pos_vel_Sat_ecef(6);

  float time = (float)timeWeek;
  for (int i = 0; i < no_sat; ++i)
  {
    uint32_t TOW = 0; //TOW = gnss_raw_measurement(i,whateverTow is)
    uint8_t L2 = 0;
    uint16_t week_No = 0;
    uint8_t L2_Flag = 0;
    uint8_t SV_Acc = 0;
    uint8_t SV_Hlth = 0;
    double T_GD = 0;
    uint16_t IODC = 0;
    double t_OC = 0;
    int8_t a_f2 = 0;
    double a_f1 = 0;
    double a_f0 = 0;
    uint8_t IODE = 0;
    double C_rs = 0;
    double delta_n = 0;
    double M_0 = 0;
    double C_uc = 0;
    double ecc = 0;
    double C_us = 0;
    double sqrt_A = 0;
    double t_OE = 0;
    double C_ic = 0;
    double Omega_0 = 0;
    double C_is = 0;
    double i_0 = 0;
    double C_rc = 0;
    double omega = 0;
    double Omega_dot = 0;
    double IDOT = 0;
    double doppler = 0;

    pos_vel_Sat_ecef = EphemerisData2Satecef(time, TOW, L2, week_No, L2_Flag,
                                             SV_Acc, SV_Hlth, T_GD, IODC, t_OC, a_f2,
                                             a_f1, a_f0, IODE, C_rs, delta_n,
                                             M_0, C_uc, ecc, C_us, sqrt_A,
                                             t_OE, C_ic, Omega_0, C_is,
                                             i_0, C_rc, omega, Omega_dot, IDOT);

    //    pseudorange, peudorange rate
    double pseudorange;
    double lambda = 2*c / (1575.42e6); // L1
    //double lambda = 2*c / (1227.60e6);            // L2
    double pseudorange_rate = lambda * doppler; //
    // record data in the row
    gnss_measurement(i, 0) = pseudorange; //gnss_raw_measurement(i,1);
    gnss_measurement(i, 1) = 0;           //-2000;
    gnss_measurement(i, 2) = pos_vel_Sat_ecef[0];
    gnss_measurement(i, 3) = pos_vel_Sat_ecef[1];
    gnss_measurement(i, 4) = pos_vel_Sat_ecef[2];
    gnss_measurement(i, 5) = pos_vel_Sat_ecef[3];
    gnss_measurement(i, 6) = pos_vel_Sat_ecef[4];
    gnss_measurement(i, 7) = pos_vel_Sat_ecef[5];
  }

  return gnss_measurement;
}

// Measurement Update
void uNavINS::MeasUpdate17(Matrix<double, Dynamic, 8> gnss_measurement)
{
  // 
  cout << "this is EKF17 filter:" << endl; 
  // determine no of satellite
  int no_sat = gnss_measurement.rows();
  MatrixXf gnss_measurement_f(no_sat,8);
  gnss_measurement_f.setZero();
  gnss_measurement_f =gnss_measurement.cast<float>();
  // Range Estimate
  MatrixXf rangeEst_m_(no_sat, 1);
  // Range rate Estimate
  MatrixXf rangeRateEst_mps_(no_sat, 1);
  // Range Error
  MatrixXf rangeErr_m(no_sat, 1);
  // Range rate Error
  MatrixXf rangeRateErr_mps(no_sat, 1);
  // Allocate measurement matrix
  MatrixXf H_(2 * no_sat, 17);
  H_.setZero();
  // Allocate measurement noise matrix
  MatrixXf R_(2 * no_sat, 2 * no_sat);
  R_.setIdentity();
  // Allocate space for Rew and Rns 
  double Rew, Rns;
  EarthRad(pEst_D_rrm_(0), &Rew, &Rns);
  // Estimate Position and Velocity in E frame 
  pEst_E_m_ = D2E(pEst_D_rrm_).cast<float>();
  vEst_E_mps_ = L2E(vEst_L_mps_, (pEst_D_rrm_).cast<float>());
  // Earth rotation vector and matrix 
  omega_ie(0) = 0.0;
  omega_ie(1) = 0.0;
  omega_ie(2) = OMEGA_DOT_EARTH;
  Omega_ie = Skew(omega_ie);  

  // loop measurement
  for (int i = 0; i < no_sat; ++i)
  {
    // predict approx range
    Vector3f delta_r;
    Vector3f pMeas_E_m_(gnss_measurement_f(i, 2), gnss_measurement_f(i, 3), gnss_measurement_f(i, 4));
    Vector3f vMeas_E_mps_(gnss_measurement_f(i, 5), gnss_measurement_f(i, 6), gnss_measurement_f(i, 7));

    delta_r = pMeas_E_m_ - pEst_E_m_;
    float approx_range =  delta_r.norm();

    // Calculate frame rotation during signal transit time using (Grove2nd:8.36)
    Matrix<double, 3, 3> T_E2I; 
    T_E2I.setZero();
    T_E2I(0, 0) = 1;
    T_E2I(0, 1) =  OMEGA_DOT_EARTH * approx_range / c;
    T_E2I(0, 2) = 0;
    T_E2I(1, 0) = -OMEGA_DOT_EARTH * approx_range / c;
    T_E2I(1, 1) = 1;
    T_E2I(1, 2) = 0;
    T_E2I(2, 0) = 0;
    T_E2I(2, 1) = 0;
    T_E2I(2, 2) = 1;

    // Predict pseudo-range using (Grove2nd:9.165)
    delta_r = T_E2I.cast<float>() * pMeas_E_m_ - pEst_E_m_;
    float range =  delta_r.norm();
    rangeEst_m_(i, 0) = range + clockBias_m_;

    // Formulate range innovation
    rangeErr_m(i, 0) = gnss_measurement_f(i, 0) - rangeEst_m_(i, 0);
    
    // Predict line of sight in ECI frame
    Matrix<float, 3, 1> u_as_E = delta_r / range;
    Matrix<double, 3, 1> u_as_L = E2L(u_as_E.cast<double>(), pEst_D_rrm_);
    // Predict pseudo-range rate using (9.165)
    float range_rate = u_as_E.transpose() * (T_E2I.cast<float>() * (vMeas_E_mps_ + Omega_ie.cast<float>() * pMeas_E_m_) 
                       - (vEst_E_mps_ + Omega_ie.cast<float>() * pEst_E_m_));
    rangeRateEst_mps_(i, 0) = range_rate  + clockRateBias_mps_;

    // Formulate range-rate innovation
    rangeRateErr_mps(i, 0) = gnss_measurement_f(i, 1) - rangeRateEst_mps_(i, 0);

    // set H in E frame 
    // H_(i, 6) = u_as_E(0);
    // H_(i, 7) = u_as_E(1);
    // H_(i, 8) = u_as_E(2);
    // H_(i, 15) = 1;
    // H_(i + no_sat, 3) = u_as_E(0);
    // H_(i + no_sat, 4) = u_as_E(1);
    // H_(i + no_sat, 5) = u_as_E(2);
    // H_(i + no_sat, 16) = 1;

    // set H in L frame - Eq.(14.128) 
    H_(i, 6) = ((float)Rns+ (float)pEst_D_rrm_(2))*u_as_L(0);
    H_(i, 7) = ((float)Rew+ (float)pEst_D_rrm_(2))*cos((float)pEst_D_rrm_(0))*u_as_L(1);
    H_(i, 8) = - (float)u_as_L(2);
    H_(i, 15) = 1;
    H_(i + no_sat, 3) = (float)u_as_L(0);
    H_(i + no_sat, 4) = (float)u_as_L(1);
    H_(i + no_sat, 5) = (float)u_as_L(2);
    H_(i + no_sat, 16) = 1;
    // set R
    R_(i, i) = pseudorangeNoiseSigma_m * pseudorangeNoiseSigma_m;
    R_(i + no_sat, i + no_sat) = pseudorangeRateNoiseSigma_mps * pseudorangeRateNoiseSigma_mps;
  }

  // Create measurement Y, as Error between Measures and Outputs
  VectorXf y(2 * no_sat, 1);
  y.setZero();
  y.segment(0, no_sat) = rangeErr_m;
  y.segment(no_sat, no_sat) = rangeRateErr_mps;

  // Innovation covariance
  MatrixXf S_(2 * no_sat, 2 * no_sat);
  S_.setZero();
  S_ = H_ * P_ * H_.transpose() + R_;

  // Kalman gain
  MatrixXf K(17, 2 * no_sat);
  K.setZero();
  K = P_ * H_.transpose() * S_.inverse();

  // Covariance update, P = (I + K * H) * P * (I + K * H)' + K * R * K'
  MatrixXf I_KH(17, 17); 
  I_KH.setZero();
  I_KH = I17 - K * H_; // temp
  P_ = I_KH * P_ * I_KH.transpose() + K * R_ * K.transpose();

  // State update, x = K * y
  VectorXf x(17, 1);
  x.setZero();
  x = K * y;

  // Pull apart x terms to update the Position, velocity, orientation, and sensor biases
  Vector3f pDeltaEst_D = x.segment(0, 3); // Position Deltas in LLA
  Vector3f vDeltaEst_L = x.segment(3, 3); // Velocity Deltas in NED
  Vector3f quatDelta = x.segment(6, 3);   // Quaternion Delta
  Vector3f aBiasDelta = x.segment(9, 3);  // Accel Bias Deltas
  Vector3f wBiasDelta = x.segment(12, 3); // Rotation Rate Bias Deltas
  float clockBiasDelta = x(15,0); // Clock Offset Deltas
  float clockBiasRateDelta = x(16,0);  // Clock Phase Deltas
  
  // Position update
  pEst_D_rrm_(2) += -pDeltaEst_D(2);
  pEst_D_rrm_(0) += pDeltaEst_D(0) / (Rew + pEst_D_rrm_(2));
  pEst_D_rrm_(1) += pDeltaEst_D(1) / (Rns + pEst_D_rrm_(2)) / cos(pEst_D_rrm_(0));

  // Velocity update
  vEst_L_mps_ += vDeltaEst_L;

  // Attitude correction
  Quaternionf dQuat_BL = Quaternionf(1.0, quatDelta(0), quatDelta(1), quatDelta(2));
  quat_BL_ = (quat_BL_ * dQuat_BL).normalized();

  // Update biases from states
  aBias_mps2_ += aBiasDelta;
  wBias_rps_ += wBiasDelta;
  clockBias_m_ += clockBiasDelta;           // clock bias
  clockRateBias_mps_ += clockBiasRateDelta; // clock rate bias
}

// Measurement Update
void uNavINS::MeasSeqUpdate17(Matrix<double, Dynamic, 8> gnss_measurement)
{
  // 
  cout << "this is EKF17 filter:" << endl; 
  // determine no of satellite
  int no_sat = gnss_measurement.rows();
  MatrixXf gnss_measurement_f(no_sat,8);
  gnss_measurement_f.setZero();
  gnss_measurement_f =gnss_measurement.cast<float>();
  // Range Estimate
  float rangeEst_m_;
  // Range rate Estimate
  float rangeRateEst_mps_;
  // Range Error
  float rangeErr_m;
  // Range rate Error
  float rangeRateErr_mps;
  // Allocate measurement matrix
  MatrixXf H_(2, 17);
  H_.setZero();
  // Allocate measurement noise matrix
  MatrixXf R_(2, 2);
  R_.setIdentity();
  // set R
  R_(0, 0) = pseudorangeNoiseSigma_m * pseudorangeNoiseSigma_m;
  R_(1, 1) = pseudorangeRateNoiseSigma_mps * pseudorangeRateNoiseSigma_mps;
  // Allocate space for Rew and Rns 
  double Rew, Rns;
  EarthRad(pEst_D_rrm_(0), &Rew, &Rns);
  // Estimate Position and Velocity in E frame 
  pEst_E_m_ = D2E(pEst_D_rrm_).cast<float>();
  vEst_E_mps_ = L2E(vEst_L_mps_, (pEst_D_rrm_).cast<float>());
  // Earth rotation vector and matrix 
  omega_ie(0) = 0.0;
  omega_ie(1) = 0.0;
  omega_ie(2) = OMEGA_DOT_EARTH;
  Omega_ie = Skew(omega_ie);  
  

  // Create measurement Y, as Error between Measures and Outputs
  VectorXf y(2, 1);
  y.setZero();

  // Innovation covariance
  MatrixXf S_(2 , 2);
  S_.setZero();

  // Kalman gain
  MatrixXf K(17, 2);
  K.setZero();

  MatrixXf I_KH(17, 17); 
  I_KH.setZero();
  

  // State update, x = K * y
  VectorXf x(17, 1);
  x.setZero();

  // loop measurement
  for (int i = 0; i < no_sat; ++i)
  {
    // predict approx range
    Vector3f delta_r;
    Vector3f pMeas_E_m_(gnss_measurement_f(i, 2), gnss_measurement_f(i, 3), gnss_measurement_f(i, 4));
    Vector3f vMeas_E_mps_(gnss_measurement_f(i, 5), gnss_measurement_f(i, 6), gnss_measurement_f(i, 7));

    delta_r = pMeas_E_m_ - pEst_E_m_;
    float approx_range =  delta_r.norm();

    // Calculate frame rotation during signal transit time using (Grove2nd:8.36)
    Matrix<double, 3, 3> T_E2I; 
    T_E2I.setZero();
    T_E2I(0, 0) = 1;
    T_E2I(0, 1) =  OMEGA_DOT_EARTH * approx_range / c;
    T_E2I(0, 2) = 0;
    T_E2I(1, 0) = -OMEGA_DOT_EARTH * approx_range / c;
    T_E2I(1, 1) = 1;
    T_E2I(1, 2) = 0;
    T_E2I(2, 0) = 0;
    T_E2I(2, 1) = 0;
    T_E2I(2, 2) = 1;

    // Predict pseudo-range using (Grove2nd:9.165)
    delta_r = T_E2I.cast<float>() * pMeas_E_m_ - pEst_E_m_;
    float range =  delta_r.norm();
    rangeEst_m_ = range + clockBias_m_;

    // Formulate range innovation
    rangeErr_m = gnss_measurement_f(i, 0) - rangeEst_m_;
    
    // Predict line of sight in ECI frame
    Matrix<float, 3, 1> u_as_E = delta_r / range;
    Matrix<double, 3, 1> u_as_L = E2L(u_as_E.cast<double>(), pEst_D_rrm_);
    // Predict pseudo-range rate using (9.165)
    float range_rate = u_as_E.transpose() * (T_E2I.cast<float>() * (vMeas_E_mps_ + Omega_ie.cast<float>() * pMeas_E_m_) 
                       - (vEst_E_mps_ + Omega_ie.cast<float>() * pEst_E_m_));
    rangeRateEst_mps_ = range_rate  + clockRateBias_mps_;

    // Formulate range-rate innovation
    rangeRateErr_mps = gnss_measurement_f(i, 1) - rangeRateEst_mps_;

    // set H in L frame - Eq.(14.128) 
    H_(0, 6) = ((float)Rns+ (float)pEst_D_rrm_(2))*u_as_L(0);
    H_(0, 7) = ((float)Rew+ (float)pEst_D_rrm_(2))*cos((float)pEst_D_rrm_(0))*u_as_L(1);
    H_(0, 8) = - (float)u_as_L(2);
    H_(0, 15) = 1;
    H_(1, 3) = (float)u_as_L(0);
    H_(1, 4) = (float)u_as_L(1);
    H_(1, 5) = (float)u_as_L(2);
    H_(1, 16) = 1;
    
    y(0) = rangeErr_m;
    y(1) = rangeRateErr_mps;
    S_ = H_ * P_ * H_.transpose() + R_;
    K = P_ * H_.transpose() * S_.inverse();
    // Covariance update, P = (I + K * H) * P * (I + K * H)' + K * R * K'
    I_KH = I17 - K * H_; // temp
    P_ = I_KH * P_ * I_KH.transpose() + K * R_ * K.transpose();
    x = K * y;
  }



  // Pull apart x terms to update the Position, velocity, orientation, and sensor biases
  Vector3f pDeltaEst_D = x.segment(0, 3); // Position Deltas in LLA
  Vector3f vDeltaEst_L = x.segment(3, 3); // Velocity Deltas in NED
  Vector3f quatDelta = x.segment(6, 3);   // Quaternion Delta
  Vector3f aBiasDelta = x.segment(9, 3);  // Accel Bias Deltas
  Vector3f wBiasDelta = x.segment(12, 3); // Rotation Rate Bias Deltas
  float clockBiasDelta = x(15,0); // Clock Offset Deltas
  float clockBiasRateDelta = x(16,0);  // Clock Phase Deltas
  
  // Position update
  pEst_D_rrm_(2) += -pDeltaEst_D(2);
  pEst_D_rrm_(0) += pDeltaEst_D(0) / (Rew + pEst_D_rrm_(2));
  pEst_D_rrm_(1) += pDeltaEst_D(1) / (Rns + pEst_D_rrm_(2)) / cos(pEst_D_rrm_(0));

  // Velocity update
  vEst_L_mps_ += vDeltaEst_L;

  // Attitude correction
  Quaternionf dQuat_BL = Quaternionf(1.0, quatDelta(0), quatDelta(1), quatDelta(2));
  quat_BL_ = (quat_BL_ * dQuat_BL).normalized();

  // Update biases from states
  aBias_mps2_ += aBiasDelta;
  wBias_rps_ += wBiasDelta;
  clockBias_m_ += clockBiasDelta;           // clock bias
  clockRateBias_mps_ += clockBiasRateDelta; // clock rate bias
}