#include <iostream>
#include "ramC.h"
#include "./EvtWnPi2.hh"
#include "EvtGenBase/EvtPDL.hh"
using namespace std;
int main(){
 EvtPDL pdl;
    pdl.read("evt.pdl");

EvtWnPi2 test;
	double pi = 3.14;
    double XM[] = {0.13957018, 0.13957018};
    double const mMu = 3;
    double GF = 1;
    RAMBOC R(2, mMu, XM);
    double sum = 0;
    int nEv = 1e5;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = R.next();
        wt *= pow(2 * pi, 4 - 3 * 3);
        EvtVector4R k = R.getV(0);
        EvtVector4R q1 = R.getV(1);
  //      EvtVector4R q2 = R.getV(2);
        EvtVector4R p = k + q1;
	EvtVector4C alpha = test.WCurrent(k,q1);
	EvtVector4C beta = alpha.conj();
	EvtComplex gama = alpha*beta;
 //       double mtr2 = 128 * pow(GF, 2)*(p * q1)*(k * q1);
 //       sum += mtr2*wt;
        if (iEv < 3) {
            cout << " Debug print at iEv=" << iEv << " =========== " << endl;
   //         cout << " k = " << k << "; m2 = " << k.mass2() << endl;
   //         cout << " q1 = " << q1 << "; m2 = " << q1.mass2() << endl;
   //         cout << " q2 = " << q2 << "; m2 = " << q2.mass2() << endl;
    //        cout << " p = " << p << "; m = " << p.mass() << endl;
            cout << " alpha = " << alpha << "beta = " << beta << "gama="<< real(gama) << "; m = " << p.mass() << endl;
	cout << "alpha p = " << alpha*p  << endl;
        }
    };
    sum /= nEv;
    double gamma = 1. / 2 / (2 * mMu) * sum;
    double gamma0 = GF * GF * pow(mMu, 5) / (192 * pow(pi, 3));
 //   std::cout << gamma << std::endl; 
    std::cout << "Works correct" << std::endl;
	
}

