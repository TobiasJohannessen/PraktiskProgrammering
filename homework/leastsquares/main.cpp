#include"../includes/matrix.h"
#include<vector>
#include<functional>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>

using mat = pp::matrix;
using vec = pp::vector;
vec lsfit(std::vector<std::function<double(double)>> f, std::vector<double> x, std::vector<double> y, std::vector<double> dy) {

    // This function performs a least squares fit of the functions in f to the data points (x, y) with uncertainties dy.
    // It returns the coefficients of the fit as a vector.

    //First, we need to check that the sizes of x, y, and dy are consistent.
    if (x.size() != y.size() || x.size() != dy.size()) {
        throw std::invalid_argument("Input vectors x, y, and dy must have the same size.");
    }

    int n = x.size();
    int m = f.size();
    mat A(n, m); // Design matrix
    vec b(n); // Right-hand side vector
    vec weights(n); // Weights for the least squares fit
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            A(i, j) = f[j](x[i]) / dy[i]; // Fill the design matrix
        }
        b[i] = y[i]/dy[i]; // Fill the right-hand side vector
    }

    // Now we need to solve the least squares problem A c = b

    vec coefficients = pp::QR::solve(A, b); // Solve for coefficients using QR decomposition
    return coefficients; // Return the coefficients of the fit
};


int main(){

    std::ifstream inputFile("data/decay.txt");
    if (!inputFile) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }
    std::ofstream outputFile("data/decay_log.txt");
    outputFile << "# t\tlog(y)\tdlny" << std::endl; // Write header to the output file
    vec x,y, dy;
    double xi, yi, dyi;
    std::string header;
    std::getline(inputFile, header); // Discard the first line
    for (std::string line; std::getline(inputFile, line);) {
        
        if (sscanf(line.c_str(), "%lf %lf %lf", &xi, &yi, &dyi) == 3) {
            x.push_back(xi);
            y.push_back(yi);
            dy.push_back(dyi);
            outputFile << xi << "\t" << std::log(yi) << "\t" << dyi / yi << std::endl; // Write the transformed data to the output file
        } else {
            std::cerr << "Error parsing line: " << line << std::endl;
        }
    }
    inputFile.close();

    std::cout << "Data to fit:" << std::endl;
    x.print("t: ");
    y.print("y: ");
    dy.print("dy: ");

    std::cout<<"\n\n";


    auto fs = std::vector<std::function<double(double)>>{
	[](double z) { return 1.0; },
	[](double z) { return -z; },
	};

    vec lny = y;
    for (int i = 0; i < lny.size(); ++i) {
        if (y[i] <= 0) {
            std::cerr << "Error: y values must be positive for logarithm." << std::endl;
            return 1;
        }
        lny[i] = std::log(y[i]); // Take the natural logarithm of y
    }
    vec dlny = dy;
    for (int i = 0; i < dlny.size(); ++i) {
        if (dy[i] <= 0) {
            std::cerr << "Error: dy values must be positive for logarithm." << std::endl;
            return 1;
        }
        dlny[i] = dy[i] / y[i]; // Calculate the relative uncertainty
    }

    std::cout << "To fit the exponential decay, we take the natural logarithm of y." << std::endl;
    std::cout << "This transforms the problem into a linear fit of ln(y) = ln(a) - b * x." << std::endl;
    std::cout << "The uncertainties dlny = dy * d/dy ln(y) = dy / y." << std::endl;
    vec coeffs = lsfit(fs, x.data, lny.data, dlny.data);
    coeffs.print("\nCoefficients:\n ");
    
    // Now we can use the coefficients to evaluate the fit at any point
    for (double z : x.data) {
        double fit_value = coeffs[0] - coeffs[1] * z; // Linear fit
        std::cout << "Fit at " << z << ": " << fit_value << std::endl;
    }


    std::ofstream fitFile("data/decay_log_fit.txt");
    std::ofstream fitFile2("data/decay_fit.txt");
    if (!fitFile) {
        std::cerr << "Error opening fit file." << std::endl;
        return 1;
    }
    fitFile << "# t\tfit_ln(y)" << std::endl;
    fitFile2 << "# t\tfit(y)" << std::endl;
    for (double x = 0; x <= 15; x += 0.01) { // Example range for fitting
        double fit_ln_y = coeffs[0] - coeffs[1] * x; // Evaluate the fit
        fitFile << x << "\t" << fit_ln_y << std::endl; // Write to file
        double fit_y = std::exp(fit_ln_y); // Convert back to y
        fitFile2 << x << "\t" << fit_y << std::endl; // Write to file
    }
    fitFile.close();
    fitFile2.close();

    std::ofstream coeffFile("data/decay_log_coefficients.txt");
    std::ofstream coeffFile2("data/decay_coefficients.txt");
    if (!coeffFile) {
        std::cerr << "Error opening coefficients file." << std::endl;
        return 1;
    }
    if (!coeffFile2) {
        std::cerr << "Error opening coefficients file." << std::endl;
        return 1;
    }
    coeffFile << coeffs[0] << "\t" << coeffs[1] << std::endl; // Write coefficients to file
    coeffFile2 << std::exp(coeffs[0]) << "\t" << coeffs[1] << std::endl; // Write coefficients to file
    coeffFile.close();
    coeffFile2.close();

    double half_life = std::log(2) / coeffs[1]; // Calculate half-life
    std::cout << "To determine half-life, calculate t_{1/2} = ln(2)/lambda: " << half_life << std::endl;
    std::cout << "Half-life of ThX, determined from fit: \n" << half_life <<  "days" <<std::endl;
    std::cout << "For reference, the half-life of ThX (Ra-224) is 3.6316(14) days." << std::endl;
    

    
    return 0;
}

