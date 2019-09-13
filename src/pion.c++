#include <iostream>
#include <fstream>
#include "ramC.h"
#include "./EvtWnPi2.hh"
#include "EvtGenBase/EvtPDL.hh"
using namespace std;


double mpi;

double spec_func_2pi(double q){
     
    EvtWnPi2 cur;
    double pi = 3.14;
    double XM[] = {mpi, mpi};
    double GF = 1;
    RAMBOC R(2, q, XM);
    double sum = 0;
    int nEv = 1e5;
    
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = R.next();
        wt *= pow(2 * pi, 4 - 2 * 3);
        EvtVector4R q1 = R.getV(0);
        EvtVector4R q2 = R.getV(1);
 //       EvtVector4R q3 = R.getV(2);

	EvtVector4C alpha = cur.WCurrent(q1,q2);
	EvtVector4C beta = alpha.conj();
	EvtComplex gama = alpha*beta/(-3 * q * q);
        sum += real(gama)*wt;
   //     double mtr2 = 128 * pow(GF, 2)*(p * q1)*(k * q1);
   //     sum += mtr2*wt;
        if (iEv < 0) {
            cout << " Debug print at iEv=" << iEv << " =========== " << endl;
           cout << " q1 = " << q1 << "; m2 = " << q1.mass2() << endl;
            cout << " q2 = " << q2 << "; m2 = " << q2.mass2() << endl;
 //           cout << " q3 = " << q3 << "; m2 = " << q3.mass2() << endl;
 //           cout << " q = " << q << "; m = " << p.mass() << endl;
            cout << "gama=" << gama << std::endl; 
        }
    };
    
        sum /= nEv;
        return sum;
    }


double spec_func_3pi(double q){
     
    EvtWnPi2 cur;
    double pi = 3.14;
    double XM[] = {mpi, mpi, mpi};
    double GF = 1;
    RAMBOC R(3, q, XM);
    double sum = 0;
    int nEv = 1e5;
    
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = R.next();
        wt *= pow(2 * pi, 4 - 3 * 3);
        EvtVector4R q1 = R.getV(0);
        EvtVector4R q2 = R.getV(1);
        EvtVector4R q3 = R.getV(2);

	EvtVector4C alpha = cur.WCurrent(q1,q2, q3);
	EvtVector4C beta = alpha.conj();
	EvtComplex gama = alpha*beta/(-3 * q * q);
        sum += real(gama)*wt;
   //     double mtr2 = 128 * pow(GF, 2)*(p * q1)*(k * q1);
   //     sum += mtr2*wt;
        if (iEv < 0) {
            cout << " Debug print at iEv=" << iEv << " =========== " << endl;
           cout << " q1 = " << q1 << "; m2 = " << q1.mass2() << endl;
            cout << " q2 = " << q2 << "; m2 = " << q2.mass2() << endl;
 //           cout << " q3 = " << q3 << "; m2 = " << q3.mass2() << endl;
 //           cout << " q = " << q << "; m = " << p.mass() << endl;
            cout << "gama=" << gama << std::endl; 
        }
    };
    
        sum /= nEv;
        return sum;
    }

double spec_func_4pi(double q){
     
    EvtWnPi2 cur;
    double pi = 3.14;
    double XM[] = {mpi, mpi, mpi, mpi};
    double GF = 1;
    RAMBOC R(4, q, XM);
    double sum = 0;
    int nEv = 1e5;
    
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = R.next();
        wt *= pow(2 * pi, 4 - 4 * 3);
        EvtVector4R q1 = R.getV(0);
        EvtVector4R q2 = R.getV(1);
        EvtVector4R q3 = R.getV(2);
        EvtVector4R q4 = R.getV(3);

	EvtVector4C alpha = cur.WCurrent(q1,q2, q3, q4);
	EvtVector4C beta = alpha.conj();
	EvtComplex gama = alpha*beta/(-3 * q * q);
        sum += real(gama)*wt;
   //     double mtr2 = 128 * pow(GF, 2)*(p * q1)*(k * q1);
   //     sum += mtr2*wt;
        if (iEv < 0) {
            cout << " Debug print at iEv=" << iEv << " =========== " << endl;
           cout << " q1 = " << q1 << "; m2 = " << q1.mass2() << endl;
            cout << " q2 = " << q2 << "; m2 = " << q2.mass2() << endl;
 //           cout << " q3 = " << q3 << "; m2 = " << q3.mass2() << endl;
 //           cout << " q = " << q << "; m = " << p.mass() << endl;
            cout << "gama=" << gama << std::endl; 
        }
    };
    
        sum /= nEv;
        return sum;
    }


double spec_func_5pi(double q){
     
    EvtWnPi2 cur;
    double pi = 3.14;
    double XM[] = {mpi, mpi, mpi, mpi, mpi};
    double GF = 1;
    RAMBOC R(5, q, XM);
    double sum = 0;
    int nEv = 1e5;
    
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = R.next();
        wt *= pow(2 * pi, 4 - 5 * 3);
        EvtVector4R q1 = R.getV(0);
        EvtVector4R q2 = R.getV(1);
        EvtVector4R q3 = R.getV(2);
        EvtVector4R q4 = R.getV(3);
        EvtVector4R q5 = R.getV(4);
        

	EvtVector4C alpha = cur.WCurrent(q1,q2, q3, q4, q5);
	EvtVector4C beta = alpha.conj();
	EvtComplex gama = alpha*beta/(-3 * q * q);
        sum += real(gama)*wt;
   //     double mtr2 = 128 * pow(GF, 2)*(p * q1)*(k * q1);
   //     sum += mtr2*wt;
        if (iEv < 0) {
            cout << " Debug print at iEv=" << iEv << " =========== " << endl;
           cout << " q1 = " << q1 << "; m2 = " << q1.mass2() << endl;
            cout << " q2 = " << q2 << "; m2 = " << q2.mass2() << endl;
 //           cout << " q3 = " << q3 << "; m2 = " << q3.mass2() << endl;
 //           cout << " q = " << q << "; m = " << p.mass() << endl;
            cout << "gama=" << gama << std::endl; 
        }
    };
    
        sum /= nEv;
        return sum;
    }


    
int main(){  
 
        EvtPDL pdl;
    pdl.read("evt.pdl");
    mpi =  EvtPDL::getMass(EvtPDL::getId("pi+"));  
    
    double qmax2 = 2;
    double qmax3 = 3;
    double qmax4 = 4;
    double qmax5 = 5;
    
    double qmin2 = 2* mpi;
    double qmin3 = 3* mpi;
    double qmin4 = 4* mpi;
    double qmin5 = 5* mpi;
    
    int qsize = 100;
    
    double qstep2 = (qmax2-qmin2)/qsize;
    double qstep3 = (qmax3-qmin3)/qsize;
    double qstep4 = (qmax4-qmin4)/qsize;
    double qstep5 = (qmax5-qmin5)/qsize;
    
    ofstream out2, out3, out4, out5; 
     out2.open("./plot2pi.txt");
     out3.open("./plot3pi.txt");
     out4.open("./plot4pi.txt");
     out5.open("./plot5pi.txt");
     
    for (int iq=1; iq<=qsize; ++iq){
        
        double qn2 = qmin2 + iq * qstep2;
        double qn3 = qmin3 + iq * qstep3;
        double qn4 = qmin4 + iq * qstep4;
        double qn5 = qmin5 + iq * qstep5;
        
        double qsum2 =  spec_func_2pi(qn2);
        double qsum3 =  spec_func_3pi(qn3);
        double qsum4 =  spec_func_4pi(qn4);
        double qsum5 =  spec_func_5pi(qn5);
        
        
        cout << iq << endl;
        out2 << qn2*qn2 << " " << qsum2  << endl;
        out3 << qn3*qn3 << " " << qsum3  << endl;
        out4 << qn4*qn4 << " " << qsum4  << endl;
        out5 << qn5*qn5 << " " << qsum5  << endl;
        
        
    };
    
    
    out2.close();
    out3.close();
    out4.close();
    out5.close();

	
}

