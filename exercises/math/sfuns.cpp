#include<cmath>
#include "sfuns.h" // Include the header file

using namespace std;


double fgamma(double x){
    ///single precision gamma function (formula from Wikipedia)
    if(x<0){
        return _PI_/sin(_PI_*x)/fgamma(1-x); // Euler's reflection formula
    };
    if(x<9){
        return fgamma(x+1)/x; // Recurrence relation;
    };
    double lnfgamma=x*log(x+1/(12*x-1/x/10))-x+log(2*_PI_/x)/2;
    
    return exp(lnfgamma);
};