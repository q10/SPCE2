#include "common.h"

const long MBIG = 1000000000;
const long MSEED = time(NULL);
const long MZ = 0;
const double FAC = (1.0 / MBIG);
double idum;

double RAN3() {
    static int inext, inextp;
    static double ma[56];
    static int iff = 0;
    long mj, mk;
    int i, ii, k;
    if (idum < 0 || iff == 0) {
        iff = 1;
        mj = MSEED - (long) (idum < 0 ? -idum : idum);
        mj = mj % MBIG;
        ma[55] = mj;
        mk = 1;
        for (i = 1; i <= 54; i++) {
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ) {
                mk += MBIG;
            }
            mj = (long) ma[ii];
        }
        for (k = 1; k <= 4; k++)
            for (i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ) {
                    ma[i] += MBIG;
                }
            }
        inext = 0;
        inextp = 31;
        idum = 1;
    }
    if (++inext == 56) inext = 1;
    if (++inextp == 56) inextp = 1;
    mj = (long) (ma[inext] - ma[inextp]);
    if (mj < MZ) {
        mj += MBIG;
    }
    ma[inext] = mj;
    return mj*FAC;

}

int RANDINT(int low, int high) {
    return (int) (RAN3()*((double) high - low)) + low;
}

double MIN(double a, double b) {
    return (a < b ? a : b);
}

double ROUND(double d) {
    return floor(d + 0.5);
}

/* The following two functions draw random numbers from a 
 * Gaussian distribution with a set mean and standard deviation.
 * Uses the Box-Mueller algorithm.
 */

double RANDGAUSS() {
    return RANDGAUSS(0.0, 1.0);
}

double RANDGAUSS(double mean, double stdev) {
    static int turn = 0;
    static double x1, x2, r;

    if (turn == 1) {
        turn = 0;
        return mean + x2 * r*stdev;
    } else {
        r = 1.0;
        while (r >= 1.0) {
            x1 = (2.0 * RAN3()) - 1.0;
            x2 = (2.0 * RAN3()) - 1.0;
            r = x1 * x1 + x2*x2;
        }
        r = sqrt(-2.0 * log(r) / r);
        turn = 1;
        return mean + x1 * r*stdev;
    }
}

double * RANDUNITVECTOR() {
    double * vector = new double[3];
    double rx, ry, r0, rz, r2 = 2.0;
    while (r2 >= 1.0) {
        rx = 1 - 2 * RAN3();
        ry = 1 - 2 * RAN3();
        r2 = rx * rx + ry*ry;
    }

    r0 = 2 * sqrt(1 - r2);
    rx *= r0;
    ry *= r0;
    rz = 1 - 2 * r2;

    vector[0] = rx;
    vector[1] = ry;
    vector[2] = rz;

    return vector;
}

string TIMESTAMP() {
    char the_date[BUF_SIZE];
    the_date[0] = '\0';
    time_t now = time(NULL);

    if (now != -1) {
        strftime(the_date, BUF_SIZE, "%Y%m%d_%H%M", localtime(&now));
    }

    return string(the_date);
}

void ASSERT(bool expression, string error_msg) {
    if (!expression) {
        cerr << "\nFATAL: " << error_msg << endl
                << "\nPREMATURELY TERMINATING PROGRAM...\n" << endl;
        abort();
    }
    return;
}
