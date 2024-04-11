#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <complex>

using namespace std;

// Calculating first derivative with the Fourier spectral method
// Example: u=exp(sin(x))
// Solution: du_dx = exp(sin(x)) * cos(x)


// Function to compute 1D FFT
fftw_complex* fft(const fftw_complex* in, const int N) {
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan forward_plan;

    // Create FFTW plan for forward transform
    forward_plan = fftw_plan_dft_1d(N, const_cast<fftw_complex*>(in), out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Execute the forward transform
    fftw_execute(forward_plan);

    // Clean up
    fftw_destroy_plan(forward_plan);

    cout << "Applying FFT" << endl;

    return out;
}

// Function to compute 1D IFFT
fftw_complex* ifft(const fftw_complex* in, const int N) {
    
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan inverse_plan;

    // Create FFTW plan for inverse transform
    inverse_plan = fftw_plan_dft_1d(N, const_cast<fftw_complex*>(in), out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Execute the inverse transform
    fftw_execute(inverse_plan);

    // Normalize the result
    for (int i = 0; i < N; ++i) {
        out[i][0] /= N;  // Normalize the real part
        out[i][1] /= N;  // Normalize the imaginary part
    }

    // Clean up
    fftw_destroy_plan(inverse_plan);

    cout << "Applying IFFT" << endl;

    return out;
}

double* fftfreq(const int size, const double sample_rate) {
    double pi = M_PI;
    // Calculate frequencies
    double* frequencies = new double[size];
    for (int i = 0; i < size; ++i) {
        if (i < size / 2) {
            frequencies[i] = 2*pi*static_cast<double>(i) / (sample_rate * size);
        }
        else {
            frequencies[i] = 2*pi*static_cast<double>(i - size) / (sample_rate * size);
        }
    }

    return frequencies;
}

int main() {
    const int Nx = 64;
    // const double dx = 1.0;
    const double pi = M_PI;
    const double L = 2*pi;
    const double dx = L/Nx;
    fftw_complex I[1]; // imaginary unit
    I[0][0] = 0;
    I[0][1] = 1;

    float x[Nx];  // x - space, 
    for (int i = 0; i < Nx; i++){
        x[i] = -L/2 + i*dx;
    }

    // u function u(x)=exp(sin(x))
    fftw_complex u[Nx];
    for (int i = 0; i < Nx; ++i) {
        u[i][0] = exp(sin(x[i])); // Real part
        u[i][1] = 0.0;   // Imaginary part
    }

    // Exact solution of u
    float u_exact[Nx]; // u_exact - exact solution of u'
    for (int i = 0; i < Nx; i++){
        u_exact[i] = exp(sin(x[i])) * cos(x[i]);
    }

    // Fourier frequencies
    double* k = fftfreq(Nx, dx);

    // Apply FFT
    fftw_complex* u_fft = fft(u, Nx);
    // cout << "\nFFT Result:" << endl;
    // for (int i = 0; i < Nx; ++i) {
    //     cout << u_fft[i][0] << " + " << u_fft[i][1] << "i" << endl;
    // }

    // solution in Fourier space
    fftw_complex u_num_fft[Nx];
    for (int i = 0; i<Nx; i++){
        u_num_fft[i][0] = k[i] * (u_fft[i][0] * I[0][0] - u_fft[i][1] * I[0][1] );
        u_num_fft[i][1] = k[i] * (u_fft[i][0] * I[0][1] - u_fft[i][1] * I[0][0] );
    }

    // real solution
    fftw_complex* u_num = ifft(u_num_fft, Nx);
    double u_num_real[Nx];
    for (int i = 0; i<Nx; i++) {
        u_num_real[i] = u_num[i][0];
    }

    // Calculate the error
    double error_fsm = 0.0;

    // Calculate the infinity norm (maximum absolute difference)
    for (int i = 0; i < Nx; ++i) {
        double absoluteDifference = std::abs(u_num_real[i] - u_exact[i]);
        error_fsm = std::max(error_fsm, absoluteDifference);
    }

    cout << "\nError = " << error_fsm << endl;

    // plot
    FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    // Set x and y labels and title
    fprintf(gnuplotPipe, "set xlabel 'x'\n");
    fprintf(gnuplotPipe, "set ylabel 'u(x)'\n");
    fprintf(gnuplotPipe, "set title '1st derivative with FSM'\n");
    // fprintf(gnuplotPipe, "set key left top\n");  // Change the legend location
    // fprintf(gnuplotPipe, "set key box\n");    // Display legend box

    // Set terminal to jpeg and specify the DPI
    fprintf(gnuplotPipe, "set term jpeg enhanced font 'arial,12' size 800,600\n");
    fprintf(gnuplotPipe, "set output 'fsm_1st_deriv.jpg'\n");

    // Plot data with different markers and set legend location
    fprintf(gnuplotPipe, "plot '-' with linespoints pt 5 title 'Analytical', '-' with lines pt 5 title 'Numerical'\n");
    
    // Data for the analytical solution
    for (size_t i = 0; i < Nx; ++i) {
        fprintf(gnuplotPipe, "%f %f\n", x[i], u_exact[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    // Data for the numerical solution
    for (size_t i = 0; i < Nx; ++i) {
        fprintf(gnuplotPipe, "%f %f\n", x[i], u_num_real[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    fflush(gnuplotPipe);

    pclose(gnuplotPipe);

    // Free memory
    fftw_free(u_fft), fftw_free(u_num);
    delete[] k, u_exact, x, u_num_real;

    return 0;
}
