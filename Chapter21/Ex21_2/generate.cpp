/*
 * High-performance Computing
 *
 * Utility program to generate a superposition of noisy sine waves.
 */
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
using namespace std;

// Pairing a frequency and an amplitude
typedef pair<double,double> T;

int main() {
    // Define the list of frequencies to use and their amplitudes.
    vector<T> freq;
    freq.push_back( T(2,  0.5) );
    freq.push_back( T(3,  0.6) );
    freq.push_back( T(13, 0.2) );
    freq.push_back( T(37, 0.4) );
    freq.push_back( T(39, 0.3) );

    double L = 10.0;
    double dx = 0.001;    // 1000Hz sampling rate
    int n = ceil(L/dx) + 1;

    // Normal distribution sampler
    const double mean = 0.0;
    const double stddev = 0.75;
    std::default_random_engine gen;
    std::normal_distribution<double> noise(mean, stddev);

    // Generate signal
    for (int i = 0; i < n; ++i) {
        double v = 0.0;
        for (auto &f : freq) {
            v += sin(f.first * 2.0 * M_PI * i*dx) * f.second + noise(gen);
        }
        cout << v << endl;
    }

    return 0;
}
