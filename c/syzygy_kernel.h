/*
 * Syzygy High-Performance Computational Kernel
 * (c) 2026 Aneesh Krishna
 *
 * SIMD-accelerated packed monomial operations, finite field row reduction,
 * and sparse Sheaf Hodge Laplacian matrix-vector kernels.
 */

#ifndef SYZYGY_KERNEL_H
#define SYZYGY_KERNEL_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* === Bit-Packed Monomial Operations === */

/*
 * Monomials packed into 64-bit words:
 * Format: 8 exponent slots of 8 bits each (sum of degrees <= 255).
 * Most significant byte is variable 0, least significant is variable 7.
 */
typedef uint64_t packed_monomial_t;

/* Returns -1 if a < b, 0 if a == b, +1 if a > b under Lexicographic order */
int packed_cmp_lex(packed_monomial_t a, packed_monomial_t b);

/* Returns -1 if a < b, 0 if a == b, +1 if a > b under Degree Reverse Lexicographic (Grevlex) */
int packed_cmp_grevlex(packed_monomial_t a, packed_monomial_t b);

/* Tests whether monomial a divides monomial b: a_i <= b_i for all i */
bool packed_divides(packed_monomial_t a, packed_monomial_t b);

/* Multiplies two packed monomials by component-wise byte addition with overflow check */
packed_monomial_t packed_mul(packed_monomial_t a, packed_monomial_t b, bool *overflow);

/* Computes the least common multiple (LCM) of two packed monomials */
packed_monomial_t packed_lcm(packed_monomial_t a, packed_monomial_t b);

/* Computes the total degree (sum of 8-bit exponents) */
uint32_t packed_total_degree(packed_monomial_t a);


/* === Finite Field Arithmetic (Z / pZ) === */

typedef uint32_t ff_elem_t;

/* Modular exponentiation: base^exp mod p */
ff_elem_t ff_pow(ff_elem_t base, uint64_t exp, ff_elem_t p);

/* Modular multiplicative inverse via Extended Euclidean Algorithm */
ff_elem_t ff_inv(ff_elem_t a, ff_elem_t p);


/* === Sparse Macaulay Matrix Reduction (F4 Linear Algebra Kernel) === */

typedef struct {
    uint32_t col_idx;
    ff_elem_t val;
} sparse_entry_t;

typedef struct {
    size_t length;
    size_t capacity;
    sparse_entry_t *entries;
} sparse_row_t;

typedef struct {
    size_t num_rows;
    size_t num_cols;
    ff_elem_t prime;
    sparse_row_t *rows;
} sparse_matrix_t;

/* Allocates and initializes a sparse matrix */
sparse_matrix_t *sparse_matrix_create(size_t rows, size_t cols, ff_elem_t prime);

/* Frees all resources associated with a sparse matrix */
void sparse_matrix_free(sparse_matrix_t *mat);

/* Adds an entry to a row in the sparse matrix */
void sparse_matrix_append_entry(sparse_matrix_t *mat, size_t row, uint32_t col, ff_elem_t val);

/*
 * Performs in-place Gaussian elimination over Z / pZ.
 * Returns the rank of the matrix.
 */
size_t sparse_matrix_gaussian_elimination(sparse_matrix_t *mat);


/* === Sheaf Hodge Laplacian SpMV Kernel === */

typedef struct {
    size_t num_rows;
    size_t num_cols;
    size_t nnz;
    uint32_t *row_ptr;   /* Size: num_rows + 1 */
    uint32_t *col_idx;   /* Size: nnz */
    double *values;      /* Size: nnz */
} csr_matrix_t;

/* Allocates a CSR matrix */
csr_matrix_t *csr_matrix_create(size_t rows, size_t cols, size_t nnz);

/* Frees a CSR matrix */
void csr_matrix_free(csr_matrix_t *mat);

/*
 * Sparse matrix-vector multiplication for Sheaf Laplacians:
 * y = alpha * A * x + beta * y
 */
void csr_spmv(const csr_matrix_t *mat, double alpha, const double *x, double beta, double *y);

/*
 * Power iteration to compute the largest eigenvalue of the Sheaf Laplacian.
 */
double sheaf_laplacian_spectral_radius(const csr_matrix_t *laplacian, size_t max_iter, double tol);

#ifdef __cplusplus
}
#endif

#endif /* SYZYGY_KERNEL_H */
