#include "ramC.h"

#include <iostream>

using namespace std;

extern "C" {
    void init_(void);
    double ram2_(double &ecm, double *xm, double *p1, double *p2);
    double ram3_(double &ecm, double *xm, double *p1, double *p2, double *p3);
    double ram4_(double &ecm, double *xm, double *p1, double *p2, double *p3, double *p4);
    double ram5_(double &ecm, double *xm, double *p1, double *p2, double *p3, double *p4, double *p5);
    double rndm2_(double);
}

TRandom r;

double rndm2_(double dummy) {
    return r.Rndm();
}

void RAMBOC::setSeed(int seed) {
    r.SetSeed(seed);
}

RAMBOC::RAMBOC(int n, const double ecm, double *xm) {
    N = n;
    ECM = ecm;
    for (int i = 0; i < N; i++) XM[i] = xm[i];
};

RAMBOC::RAMBOC(double ecm, double m1, double m2) {
    N = 2;
    ECM = ecm;
    XM[0] = m1;
    XM[1] = m2;
}

RAMBOC::RAMBOC(double ecm, double m1, double m2, double m3) {
    N = 3;
    ECM = ecm;
    XM[0] = m1;
    XM[1] = m2;
    XM[2] = m3;
}

RAMBOC::RAMBOC(double ecm, double m1, double m2, double m3, double m4) {
    N = 4;
    ECM = ecm;
    XM[0] = m1;
    XM[1] = m2;
    XM[2] = m3;
    XM[3] = m4;
}

RAMBOC::RAMBOC(double ecm, double m1, double m2, double m3, double m4, double m5) {
    N = 5;
    ECM = ecm;
    XM[0] = m1;
    XM[1] = m2;
    XM[2] = m3;
    XM[3] = m4;
    XM[4] = m5;
}

double RAMBOC::next() {
    if (N == 2) {
        WT = ram2_(ECM, XM, p1, p2);
        P1.set(p1[3], p1[0], p1[1], p1[2]);
        P2.set(p2[3], p2[0], p2[1], p2[2]);
    } else if (N == 3) {
        WT = ram3_(ECM, XM, p1, p2, p3);
        P1.set(p1[3], p1[0], p1[1], p1[2]);
        P2.set(p2[3], p2[0], p2[1], p2[2]);
        P3.set(p3[3], p3[0], p3[1], p3[2]);
    } else if (N == 4) {
        WT = ram4_(ECM, XM, p1, p2, p3, p4);
        P1.set(p1[3], p1[0], p1[1], p1[2]);
        P2.set(p2[3], p2[0], p2[1], p2[2]);
        P3.set(p3[3], p3[0], p3[1], p3[2]);
        P4.set(p4[3], p4[0], p4[1], p4[2]);
    } else if (N == 5) {
        WT = ram5_(ECM, XM, p1, p2, p3, p4, p5);
        P1.set(p1[3], p1[0], p1[1], p1[2]);
        P2.set(p2[3], p2[0], p2[1], p2[2]);
        P3.set(p3[3], p3[0], p3[1], p3[2]);
        P4.set(p4[3], p4[0], p4[1], p4[2]);
        P5.set(p5[3], p5[0], p5[1], p5[2]);
    };
    return WT;
};

EvtVector4R RAMBOC::getV(int i) {
    if (i == 0) return P1;
    else if (i == 1) return P2;
    else if (i == 2) return P3;
    else if (i == 3) return P4;
    else if (i == 4) return P5;
    else {
        cout << "ramC: not supported yet" << endl;
        abort();
    }
};

