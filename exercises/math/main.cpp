#include <iostream>
#include <string>
#include <cmath>
#include "sfuns.cpp"

using namespace std;






int main(){

	cout << "Simple math operations:" << endl;

	
	double sqrt2 = sqrt(2.0);

	cout << "The square root of 2 equals " << sqrt2 << ". Should be 1.4142145 "<< endl;

	double root5_2 = pow(2, 1/5.0);

	cout << "The fifth root of 2 equals " << root5_2 << ". Should be 1.1486983 " << endl;

	double e_pi = pow(_E_,_PI_);

	cout << "e to the power of pi equals " << e_pi << ". Should be 23.140692 "<< endl;

	double pi_e = pow(_PI_,_E_);

	cout << "pi to the power of pi equals " << pi_e << ". Should be 22.459157 "<< endl;

	cout << endl;

	cout << "Gamma functions of 1-10:" << endl;

	for(int i=1; i<=10; i++){
		cout << "Γ(" << i << ") = " << fgamma(i) << endl;
	};


	return 0;
};
