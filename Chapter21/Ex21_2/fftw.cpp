/*
 * High-performance Computing
 *
 * Decompose noisy signal to generate spectrum.
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

#include <fftw3.h>

int main() {
    std::vector<double> s;
    double val;

    // Read signal into s until end of file
    ifstream sig("signal.txt");
    while (!sig.eof()) {
        sig >> val;
        s.push_back(val);
    }
    sig.close();

    const int n = s.size();
    const int n_out = n / 2 + 1;
    const double srate = 1000.0;        // Sample rate = 1000Hz

    // Set up output array of (fftw) complex numbers
    fftw_complex* tmp = new fftw_complex[n_out];

    // Plan the transform. We have real 1D data, so we use the real-to-complex
    // (r2c) transform. The FFTW_ESTIMATE tells FFTW how much effort to put into
    // creating an efficient plan. Since we will only do the transform once, we
    // just ask for an estimate.
    fftw_plan p = fftw_plan_dft_r2c_1d(s.size(), &s[0], tmp, FFTW_ESTIMATE);

    // Execute the plan. This actually performs the transform.
    fftw_execute(p);

    // Clean up
    fftw_destroy_plan(p);

    // The size of each frequency bin is sample_rate / num_samples
    double fbinsize = srate/n;

    // Output the power = 2*|w|, where w is a complex frequency
    std::vector<double> out(n_out, 0.0);
    for (int i = 0; i < n_out; ++i) {
        out[i] = sqrt((tmp[i][0] * tmp[i][0] + tmp[i][1] * tmp[i][1]))/(n);
        if (i > 0) {
            out[i] *= 2.0;
        }
        cout << i*fbinsize << "  " << out[i] << endl;
    }

    return 0;
}
