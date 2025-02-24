#include<cmath>
#include "sfuns.h" // Include the header file

using namespace std;


double fgamma(double x){
    ///single precision gamma function (formula from Wikipedia)
    if(x<0){
        return M_PI/sin(M_PI*x)/fgamma(1-x); // Euler's reflection formula
    };
    if(x<9){
        return fgamma(x+1)/x; // Recurrence relation;
    };
    double lnfgamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
    
    return exp(lnfgamma);
};


double lngamma(double x){
    ///single precision gamma function (formula from Wikipedia)
    if(x <= 0){ 
        return NAN;
    };
    if(x < 9){ 
        return lngamma(x+1) - std::log(x);
    };
    double lnfgamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
    
    return lnfgamma;
};


