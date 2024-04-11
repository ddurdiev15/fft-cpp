#include <iostream>
#include <fftw3.h>

int main() {
    const int size = 8;
    const double sample_rate = 1.0; // Replace this with your actual sample rate

    double input[size] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };

    fftw_complex* output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (size));

    fftw_plan forward_plan = fftw_plan_dft_r2c_1d(size, input, output, FFTW_ESTIMATE);
    fftw_plan backward_plan = fftw_plan_dft_c2r_1d(size, output, input, FFTW_ESTIMATE);

    fftw_execute(forward_plan);

    // Calculate frequencies
    double* frequencies = new double[size];
    for (int i = 0; i < size; ++i) {
        if (i < size / 2) {
            frequencies[i] = static_cast<double>(i) / (sample_rate * size);
        }
        else {
            frequencies[i] = static_cast<double>(i - size) / (sample_rate * size);
        }
    }

    // Print the frequencies and corresponding FFT result
    std::cout << "\nFrequencies:" << std::endl;
    for (int i = 0; i < size; ++i) {
        std::cout << frequencies[i] << " ";
    }
    std::cout << std::endl;

    // Print the FFT result
    std::cout << "\nFFT Result:" << std::endl;
    for (int i = 0; i < size; ++i) {
        std::cout << output[i][0] << " + " << output[i][1] << "i" << std::endl;
    }

    fftw_execute(backward_plan);

    // Normalize the result after the IFFT
    for (int i = 0; i < size; ++i) {
        input[i] /= size;
    }

    // Print the back-transformed result
    std::cout << "\nBack-transformed Result:" << std::endl;
    for (int i = 0; i < size; ++i) {
        std::cout << input[i] << " ";
    }
    std::cout << std::endl;

    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
    fftw_free(output);
    delete[] frequencies;

    return 0;
}
