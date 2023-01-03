/*
 * Benchmark and verification harness for Syzygy C acceleration kernel.
 * (c) 2026 Aneesh Krishna
 */

#include "syzygy_kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

static double get_time_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(void) {
    printf("=== Syzygy C Kernel Benchmark & Verification ===\n");

    /* 1. Verify packed monomial operations */
    printf("[1/3] Testing packed monomial SIMD operations...\n");
    packed_monomial_t m1 = 0x0302010000000000ULL; /* x_0^3 x_1^2 x_2^1 */
    packed_monomial_t m2 = 0x0201010000000000ULL; /* x_0^2 x_1^1 x_2^1 */

    assert(packed_cmp_lex(m1, m2) > 0);
    assert(packed_total_degree(m1) == 6);
    assert(packed_total_degree(m2) == 4);
    assert(packed_divides(m2, m1) == true);
    assert(packed_divides(m1, m2) == false);

    bool overflow = false;
    packed_monomial_t m_mul = packed_mul(m1, m2, &overflow);
    assert(!overflow);
    assert(packed_total_degree(m_mul) == 10);

    /* 2. Benchmark finite field sparse Gaussian elimination */
    printf("[2/3] Benchmarking sparse Macaulay matrix elimination (F4 kernel)...\n");
    size_t rows = 500;
    size_t cols = 600;
    ff_elem_t p = 65521; /* Largest 16-bit prime */
    sparse_matrix_t *mat = sparse_matrix_create(rows, cols, p);

    srand(42);
    for (size_t r = 0; r < rows; ++r) {
        /* Add ~15 non-zero entries per row */
        for (int k = 0; k < 15; ++k) {
            uint32_t col = (uint32_t)(rand() % cols);
            ff_elem_t val = (ff_elem_t)(1 + (rand() % (p - 1)));
            sparse_matrix_append_entry(mat, r, col, val);
        }
    }

    double t0 = get_time_sec();
    size_t rank = sparse_matrix_gaussian_elimination(mat);
    double t1 = get_time_sec();

    printf("      Eliminated %zu x %zu sparse matrix in %.4f seconds (rank = %zu)\n",
           rows, cols, t1 - t0, rank);
    sparse_matrix_free(mat);

    /* 3. Sheaf Laplacian spectral radius */
    printf("[3/3] Testing Sheaf Hodge Laplacian SpMV & spectral radius...\n");
    size_t n = 100;
    csr_matrix_t *lap = csr_matrix_create(n, n, n * 3);
    lap->row_ptr[0] = 0;
    size_t idx = 0;
    for (size_t i = 0; i < n; ++i) {
        lap->col_idx[idx] = (uint32_t)i;
        lap->values[idx] = 2.0;
        idx++;
        if (i > 0) {
            lap->col_idx[idx] = (uint32_t)(i - 1);
            lap->values[idx] = -1.0;
            idx++;
        }
        if (i + 1 < n) {
            lap->col_idx[idx] = (uint32_t)(i + 1);
            lap->values[idx] = -1.0;
            idx++;
        }
        lap->row_ptr[i + 1] = (uint32_t)idx;
    }
    lap->nnz = idx;

    double lambda_max = sheaf_laplacian_spectral_radius(lap, 200, 1e-6);
    printf("      Estimated Sheaf Laplacian spectral radius: %.6f\n", lambda_max);
    assert(lambda_max > 0.0 && lambda_max <= 4.01);
    csr_matrix_free(lap);

    printf("PASSED: All C acceleration kernel assertions verified.\n");
    return 0;
}
