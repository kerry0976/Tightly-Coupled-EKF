#ifndef GNSS_MEASUREMENT
#define GNSS_MEASUREMENT

#include <stdint.h>
#include <string>
#include <vector>
using namespace std;

#pragma once

struct GNSS_raw_measurement
{
    double AODO;          // the age of data offset, in seconds.
    double Cic;           // Amplitude of the Cosine Harmonic Correction Term to the Angle of Inclination
    double Cis;           // Amplitude of the Sine Harmonic Correction Term to the Angle of Inclination
    double Crc;           // Amplitude of the Cosine Harmonic Correction Term to the Orbit Radius
    double Crs;           // Amplitude of the Sine Harmonic Correction Term to the Orbit Radius
    double Cuc;           // Amplitude of the Cosine Harmonic Correction Term to the Argument of Latitude
    double Cus;           // Amplitude of the Sine Harmonic Correction Term to the Argument of Latitude
    bool FIT;             // ?
    double IDOT;          // Rate of Inclination Angle
    int IODC;             // Issue of Data, Clock (10 Bits) [Units: N/A]
    int IODE;             // Issue of Data (8 Bits) [Units: N/A]
    int L2;               // Code on L2 (2 Bits) [Units: N/A]
    int L2P;              // L2 P Data Flag (1 Bit) [Units: Discrete]
    double M0;            // Mean Anomaly at Reference Time
    double Omega0;        // Longitude of Ascending Node of Orbit Plane at Weekly Epoch
    double Omegad;        // Rate of Right Ascension
    int TOW17;            // Time-of-week
    double Tgd;           // (8 Bits / Two's Complement with sign on MSB) [Units: Seconds]
    double WN;            // GPS Week Number (10 Bits) [units: Week]
    double af0;           // (22 Bits / Two's Complement with sign on MSB) [Units: Seconds]
    double af1;           // (16 Bits / Two's Complement with sign on MSB) [Units: Sec/Sec]
    double af2;           // (8 Bits / Two's Complement with sign on MSB) [Units: Sec/Sec^2]
    string constellation; // type of constellation, e.g., GPS, GLONASS, BEIDOU
    double deltan;        // Mean Motion Difference From Computed Value
    double doppler;       // Doppler shift measurement, note: by multiplying frequency, we get psedu-range rate
    double e;             // Eccentricity
    bool frame1;          // Validity of subframe 1
    bool frame2;          // Validity of subframe 2
    bool frame3;          // Validity of subframe 3
    int gnssid;           // Satellite ID number
    int hlth;             // Satellite Vehicle Health (6 Bits) [Units: Discretes]
    double i0;            // Inclination at Reference Time
    double omega;         // Argument of Perigee
    double pseudorange;   // pseudorange (m)
    double sqrtA;         // Square Root off the Semi-Major Axis
    double timestamp;     // current seconds
    double toc;           // (16 Bits) [Units: Seconds]
    double toe;           // Reference Time Ephemeris
    int ura;              // Satellite Vehicle Accuracy (4 Bits) [Units: N/A], binary 
};

#endif