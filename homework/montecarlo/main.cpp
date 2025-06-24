#include <cmath>
#include <functional>
#include <iostream>
#include<random>
#include <vector>
#include <iomanip>
#include<fstream>

std::pair<double, double> montecarlo(std::function<double(std::vector<double>)> f, std::vector<double> a, std::vector<double> b, int N) {
    int dim = a.size();
    double volume = 1.0;
    for (int i = 0; i < dim; ++i) {
        volume *= (b[i] - a[i]);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::vector<std::uniform_real_distribution<>> distributions;
    for (int i = 0; i < dim; ++i) {
        distributions.emplace_back(a[i], b[i]);
    }

    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < N; ++i) {
        std::vector<double> point(dim);
        for (int j = 0; j < dim; ++j) {
            point[j] = distributions[j](gen);
        }

        double fx = f(point);
        sum += fx;
        sum2 += fx * fx;
    }

    double mean = sum / N;
    double variance = (sum2 / N) - (mean * mean);
    
    return {mean * volume, std::sqrt(variance / N) * volume};
}


double unit_circle(std::vector<double> p){
    double x = p[0];
    double y = p[1];
    return (x * x + y * y <= 1.0) ? 1.0 : 0.0; // Indicator function for unit circle
}


double difficult_integral(std::vector<double> p) {
    double x = p[0];
    double y = p[1];
    double z = p[2];
    double denominator = 1.0 - cos(x) * cos(y) * cos(z);
    return 1.0 / (denominator * M_PI * M_PI * M_PI);
}
int main(){

    std::vector<bool> test_results;
    bool all_tests_passed = true;

    auto run_test = [&](const std::string& name, long double result, long double expected, long double tol = 1e-4) {
        bool pass = std::fabs(result - expected) < tol;
        std::cout << std::fixed << std::setprecision(10);
        std::cout << name << ":\n"
                  << "  Result:   " << result << "\n"
                  << "  Expected: " << expected << "\n"
                  << "  " << (pass ? "PASSED" : "FAILED") << "\n\n";
        test_results.push_back(pass);
        if (!pass) all_tests_passed = false;
    };
    // Define the integration limits for a 2D unit square

    std::cout << "--------------------TESTING MONTE CARLO INTEGRATION-----------------:\n";
    std::cout << "Estimating the area of a unit circle using Monte Carlo method.\n\n";
    std::vector<double> a = {0.0, 0.0};
    std::vector<double> b = {1.0, 1.0};

    
    for (int N= 1000; N <= 100000; N *= 10) {
        // Perform Monte Carlo integration
        auto result = montecarlo(unit_circle, a, b, N);

        // Output the result
        std::cout << "N = " << N << ": Estimated integral = " << result.first 
                  << ", Estimated error = " << result.second << std::endl;
        run_test("How close?", result.first, M_PI / 4.0, 1e-2);
    }
    //Gather data points to estimate error function:
    std::ofstream results_file("data/circle_errors.txt");
    for (int N = 1000; N <= 1000000; N += N/10) {
        auto result = montecarlo(unit_circle, a, b, N);
        results_file << N << "\t" << result.first << "\t" << result.second << "\n";
        
    }

    // Test the difficult integral

    std::cout << "\n--------------------TESTING DIFFICULT INTEGRAL-----------------:\n";
    double difficult_integral_result = 1.3932039296856768591842462603255;
    std::vector<double> a_difficult = {0.0, 0.0, 0.0};
    std::vector<double> b_difficult = {M_PI, M_PI, M_PI};
    for (int N = 1000; N <= 1000000; N *= 10) {
        auto result_difficult = montecarlo(difficult_integral, a_difficult, b_difficult, N);
        std::cout << "N = " << N << ": Estimated integral = " << result_difficult.first 
                  << ", Estimated error = " << result_difficult.second << std::endl;
        run_test("Difficult integral", result_difficult.first, difficult_integral_result, 1e-2);
    }
    
    return 0;
}