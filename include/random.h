/* 
 * File:   random.h
 * Author: BENSON J MA
 *
 * Created on August 3, 2011, 4:50 PM
 */

#pragma once
#ifndef RANDOM_H
#define	RANDOM_H

#define BUF_SIZE 1024
#define D_INFINITY std::numeric_limits<double>::infinity()

using namespace std;

typedef complex<double> dcomplex;

double RAN3();
int RANDINT(int low, int high);
double MIN(double a, double b);
double ROUND(double d);
double RANDGAUSS(double mean, double stdev);
double RANDGAUSS();
double * RANDUNITVECTOR();

string TIMESTAMP();
void ASSERT(bool expression, string error_msg);
long long int timeval_diff(struct timeval *end_time, struct timeval *start_time);

template <typename T> string STRING(T tval) {
    stringstream out;
    out << tval;
    return out.str();
}

inline double DEG2RAD(double degrees) {
    return degrees * M_PI / 180.0;
}

#endif	/* RANDOM_H */
