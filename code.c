// Your optimized implementation. This is the only file you submit.
//
// Rules (submissions that break them are disqualified):
//   - Keep the fft_student signature from fft.h. Any helper functions or
//     file-scope variables you add must be declared `static`.
//   - Results must match the reference to double precision (relative error <= 1e-8).
//   - Only math, memory (memcpy/memset/malloc/free) and printf functions may be called:
//     no threads, timers, file I/O, or system calls.
//   - Do not keep results between calls; every call gets new input.
//   - No #pragma GCC optimize or __attribute__((optimize)): your code is compiled at -O0.
//   - Individual work only: do not discuss your solution, work in pairs, or share code.
//
// Test with:  make test        (or: make sizes)
#include <math.h>
#include "fft.h"

#define PI 3.14159265358979323846

void fft_student(double in[N][N],
                 double re[N][N],
                 double im[N][N],
                 double amp[N][N])
{
    double c[N][N];
    double s[N][N];
    double tmp_re[N][N];
    double tmp_im[N][N];

    const double step = 2.0 * PI / N;
    const double scale = 1.0 / N;

    int k, x, y;

    /*
     * Precompute twiddle factors.
     * c[k][x] = cos(2*pi*k*x/N)
     * s[k][x] = sin(2*pi*k*x/N)
     */
    for (k = 0; k < N; k++) {
        for (x = 0; x < N; x++) {
            int index = (k * x) % N;
            double angle = step * index;

            c[k][x] = cos(angle);
            s[k][x] = sin(angle);
        }
    }

    /*
     * First pass: transform each row.
     * Store the complex intermediate result.
     */
    for (y = 0; y < N; y++) {
        for (k = 0; k < N; k++) {
            double sum_re = 0.0;
            double sum_im = 0.0;

            for (x = 0; x < N; x++) {
                double value = in[y][x];

                sum_re += value * c[k][x];
                sum_im -= value * s[k][x];
            }

            tmp_re[y][k] = sum_re;
            tmp_im[y][k] = sum_im;
        }
    }

    /*
     * Second pass: transform columns.
     * Apply the required overall normalization of 1/N.
     */
    for (k = 0; k < N; k++) {
        for (x = 0; x < N; x++) {
            double sum_re = 0.0;
            double sum_im = 0.0;

            for (y = 0; y < N; y++) {
                double a = tmp_re[y][x];
                double b = tmp_im[y][x];
                double co = c[k][y];
                double si = s[k][y];

                sum_re += a * co + b * si;
                sum_im += b * co - a * si;
            }

            re[k][x] = sum_re * scale;
            im[k][x] = sum_im * scale;

            amp[k][x] = sqrt(re[k][x] * re[k][x]
                           + im[k][x] * im[k][x]);
        }
    }
}
