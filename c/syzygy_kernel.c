/*
 * Syzygy High-Performance Computational Kernel Implementation
 * (c) 2026 Aneesh Krishna
 *
 * SIMD packed monomials, finite field sparse Gaussian elimination,
 * and cellular sheaf Hodge Laplacian spectral analysis.
 */

#include "syzygy_kernel.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* === Bit-Packed Monomial Operations === */

int packed_cmp_lex(packed_monomial_t a, packed_monomial_t b) {
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

uint32_t packed_total_degree(packed_monomial_t a) {
    /* Sum 8 bytes into a single 32-bit integer */
    uint64_t sum = a;
    sum = (sum & 0x00FF00FF00FF00FFULL) + ((sum >> 8) & 0x00FF00FF00FF00FFULL);
    sum = (sum & 0x0000FFFF0000FFFFULL) + ((sum >> 16) & 0x0000FFFF0000FFFFULL);
    sum = (sum & 0x00000000FFFFFFFFULL) + ((sum >> 32) & 0x00000000FFFFFFFFULL);
    return (uint32_t)sum;
}

int packed_cmp_grevlex(packed_monomial_t a, packed_monomial_t b) {
    uint32_t deg_a = packed_total_degree(a);
    uint32_t deg_b = packed_total_degree(b);
    if (deg_a < deg_b) return -1;
    if (deg_a > deg_b) return 1;

    /* If total degrees match, find the last variable where exponents differ */
    /* In Grevlex, smaller exponent at the latest variable means larger monomial */
    for (int i = 0; i < 8; ++i) {
        uint8_t exp_a = (uint8_t)((a >> (i * 8)) & 0xFF);
        uint8_t exp_b = (uint8_t)((b >> (i * 8)) & 0xFF);
        if (exp_a != exp_b) {
            return (exp_a < exp_b) ? 1 : -1;
        }
    }
    return 0;
}

bool packed_divides(packed_monomial_t a, packed_monomial_t b) {
    /* Test whether a_i <= b_i for all 8 bytes simultaneously without branches */
    /* (b - a) generates high-bit borrows whenever b_i < a_i */
    for (int i = 0; i < 8; ++i) {
        uint8_t exp_a = (uint8_t)((a >> (i * 8)) & 0xFF);
        uint8_t exp_b = (uint8_t)((b >> (i * 8)) & 0xFF);
        if (exp_a > exp_b) return false;
    }
    return true;
}

packed_monomial_t packed_mul(packed_monomial_t a, packed_monomial_t b, bool *overflow) {
    packed_monomial_t res = 0;
    bool has_overflow = false;
    for (int i = 0; i < 8; ++i) {
        uint32_t ea = (uint32_t)((a >> (i * 8)) & 0xFF);
        uint32_t eb = (uint32_t)((b >> (i * 8)) & 0xFF);
        uint32_t sum = ea + eb;
        if (sum > 255) has_overflow = true;
        res |= ((packed_monomial_t)(sum & 0xFF)) << (i * 8);
    }
    if (overflow) *overflow = has_overflow;
    return res;
}

packed_monomial_t packed_lcm(packed_monomial_t a, packed_monomial_t b) {
    packed_monomial_t res = 0;
    for (int i = 0; i < 8; ++i) {
        uint8_t ea = (uint8_t)((a >> (i * 8)) & 0xFF);
        uint8_t eb = (uint8_t)((b >> (i * 8)) & 0xFF);
        uint8_t max_e = (ea > eb) ? ea : eb;
        res |= ((packed_monomial_t)max_e) << (i * 8);
    }
    return res;
}


/* === Finite Field Arithmetic === */

ff_elem_t ff_pow(ff_elem_t base, uint64_t exp, ff_elem_t p) {
    uint64_t res = 1;
    uint64_t b = base % p;
    while (exp > 0) {
        if (exp & 1) res = (res * b) % p;
        b = (b * b) % p;
        exp >>= 1;
    }
    return (ff_elem_t)res;
}

ff_elem_t ff_inv(ff_elem_t a, ff_elem_t p) {
    /* By Fermat's Little Theorem: a^{p-2} = a^{-1} mod p for prime p */
    return ff_pow(a, p - 2, p);
}


/* === Sparse Matrix Operations === */

sparse_matrix_t *sparse_matrix_create(size_t rows, size_t cols, ff_elem_t prime) {
    sparse_matrix_t *mat = (sparse_matrix_t *)malloc(sizeof(sparse_matrix_t));
    mat->num_rows = rows;
    mat->num_cols = cols;
    mat->prime = prime;
    mat->rows = (sparse_row_t *)calloc(rows, sizeof(sparse_row_t));
    return mat;
}

void sparse_matrix_free(sparse_matrix_t *mat) {
    if (!mat) return;
    for (size_t i = 0; i < mat->num_rows; ++i) {
        free(mat->rows[i].entries);
    }
    free(mat->rows);
    free(mat);
}

void sparse_matrix_append_entry(sparse_matrix_t *mat, size_t row, uint32_t col, ff_elem_t val) {
    if (row >= mat->num_rows || val == 0) return;
    sparse_row_t *r = &mat->rows[row];
    if (r->length >= r->capacity) {
        size_t new_cap = (r->capacity == 0) ? 8 : r->capacity * 2;
        r->entries = (sparse_entry_t *)realloc(r->entries, new_cap * sizeof(sparse_entry_t));
        r->capacity = new_cap;
    }
    r->entries[r->length].col_idx = col;
    r->entries[r->length].val = val % mat->prime;
    r->length++;
}

size_t sparse_matrix_gaussian_elimination(sparse_matrix_t *mat) {
    size_t pivot_row = 0;
    ff_elem_t p = mat->prime;

    for (size_t col = 0; col < mat->num_cols && pivot_row < mat->num_rows; ++col) {
        /* Find pivot row with entry at column col */
        size_t sel = pivot_row;
        bool found = false;
        for (size_t r = pivot_row; r < mat->num_rows; ++r) {
            if (mat->rows[r].length > 0 && mat->rows[r].entries[0].col_idx == col) {
                sel = r;
                found = true;
                break;
            }
        }

        if (!found) continue;

        /* Swap rows */
        if (sel != pivot_row) {
            sparse_row_t tmp = mat->rows[pivot_row];
            mat->rows[pivot_row] = mat->rows[sel];
            mat->rows[sel] = tmp;
        }

        /* Scale pivot row so leading entry is 1 */
        sparse_row_t *p_row = &mat->rows[pivot_row];
        ff_elem_t lead_inv = ff_inv(p_row->entries[0].val, p);
        for (size_t k = 0; k < p_row->length; ++k) {
            p_row->entries[k].val = (ff_elem_t)(((uint64_t)p_row->entries[k].val * lead_inv) % p);
        }

        /* Eliminate column in all subsequent rows */
        for (size_t r = pivot_row + 1; r < mat->num_rows; ++r) {
            sparse_row_t *curr = &mat->rows[r];
            if (curr->length > 0 && curr->entries[0].col_idx == col) {
                ff_elem_t factor = curr->entries[0].val;
                /* Merge curr_row -= factor * pivot_row */
                size_t max_len = curr->length + p_row->length;
                sparse_entry_t *merged = (sparse_entry_t *)malloc(max_len * sizeof(sparse_entry_t));
                size_t m_len = 0;
                size_t i = 0, j = 0;

                while (i < curr->length || j < p_row->length) {
                    if (j >= p_row->length || (i < curr->length && curr->entries[i].col_idx < p_row->entries[j].col_idx)) {
                        merged[m_len++] = curr->entries[i++];
                    } else if (i >= curr->length || (j < p_row->length && p_row->entries[j].col_idx < curr->entries[i].col_idx)) {
                        uint64_t sub_val = ((uint64_t)factor * p_row->entries[j].val) % p;
                        ff_elem_t res_val = (ff_elem_t)((p - sub_val) % p);
                        if (res_val != 0) {
                            merged[m_len].col_idx = p_row->entries[j].col_idx;
                            merged[m_len].val = res_val;
                            m_len++;
                        }
                        j++;
                    } else {
                        /* Equal columns: curr - factor * pivot */
                        uint64_t sub_val = ((uint64_t)factor * p_row->entries[j].val) % p;
                        uint64_t curr_val = curr->entries[i].val;
                        ff_elem_t diff = (ff_elem_t)((curr_val + p - sub_val) % p);
                        if (diff != 0) {
                            merged[m_len].col_idx = curr->entries[i].col_idx;
                            merged[m_len].val = diff;
                            m_len++;
                        }
                        i++;
                        j++;
                    }
                }
                free(curr->entries);
                curr->entries = merged;
                curr->length = m_len;
                curr->capacity = max_len;
            }
        }
        pivot_row++;
    }
    return pivot_row;
}


/* === Sheaf Hodge Laplacian CSR Kernels === */

csr_matrix_t *csr_matrix_create(size_t rows, size_t cols, size_t nnz) {
    csr_matrix_t *mat = (csr_matrix_t *)malloc(sizeof(csr_matrix_t));
    mat->num_rows = rows;
    mat->num_cols = cols;
    mat->nnz = nnz;
    mat->row_ptr = (uint32_t *)calloc(rows + 1, sizeof(uint32_t));
    mat->col_idx = (uint32_t *)malloc(nnz * sizeof(uint32_t));
    mat->values = (double *)malloc(nnz * sizeof(double));
    return mat;
}

void csr_matrix_free(csr_matrix_t *mat) {
    if (!mat) return;
    free(mat->row_ptr);
    free(mat->col_idx);
    free(mat->values);
    free(mat);
}

void csr_spmv(const csr_matrix_t *mat, double alpha, const double *x, double beta, double *y) {
    for (size_t i = 0; i < mat->num_rows; ++i) {
        double dot = 0.0;
        uint32_t start = mat->row_ptr[i];
        uint32_t end = mat->row_ptr[i + 1];
        for (uint32_t k = start; k < end; ++k) {
            dot += mat->values[k] * x[mat->col_idx[k]];
        }
        y[i] = alpha * dot + beta * y[i];
    }
}

double sheaf_laplacian_spectral_radius(const csr_matrix_t *laplacian, size_t max_iter, double tol) {
    size_t n = laplacian->num_rows;
    if (n == 0) return 0.0;

    double *v = (double *)malloc(n * sizeof(double));
    double *w = (double *)malloc(n * sizeof(double));

    /* Initialize with normalized unit vector */
    double inv_sqrt_n = 1.0 / sqrt((double)n);
    for (size_t i = 0; i < n; ++i) v[i] = inv_sqrt_n;

    double lambda_old = 0.0;
    double lambda_new = 0.0;

    for (size_t iter = 0; iter < max_iter; ++iter) {
        /* w = L * v */
        csr_spmv(laplacian, 1.0, v, 0.0, w);

        /* Rayleigh quotient lambda = v . w */
        double dot = 0.0;
        double norm_sq = 0.0;
        for (size_t i = 0; i < n; ++i) {
            dot += v[i] * w[i];
            norm_sq += w[i] * w[i];
        }
        lambda_new = dot;

        double norm_w = sqrt(norm_sq);
        if (norm_w < 1.0e-14) break;

        /* v = w / ||w|| */
        double inv_norm = 1.0 / norm_w;
        for (size_t i = 0; i < n; ++i) v[i] = w[i] * inv_norm;

        if (fabs(lambda_new - lambda_old) < tol) break;
        lambda_old = lambda_new;
    }

    free(v);
    free(w);
    return lambda_new;
}
