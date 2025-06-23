#include"../includes/matrix.h"
#include"../includes/ODE.h"
#include<iostream>
#include<fstream>



vec f_exp(double x, const vec& y){return 1.0/10* y;};


int main(){
    double x_init, x_final;
    x_init = 0.1; x_final = 10;
    vec y_init = vec(1);
    for (int i = 0; i < y_init.size(); ++i) {
        y_init[i] = 1.0 / (i + 1); // Initial values for y
    }

    auto [xs, ys] =  driver(f_exp, x_init, x_final, y_init, 0.05, 0.01, 0.01);

    // Print the results
    std::ofstream output_file("exponential.txt");

    for (size_t i = 0; i < xs.size(); ++i) {
        output_file << xs[i] << " ";
        for (int j = 0; j < ys[i].size(); ++j) {
            output_file << ys[i][j] << " ";
        }
        output_file << "\n";
    }
    output_file.close();
return 0;
}