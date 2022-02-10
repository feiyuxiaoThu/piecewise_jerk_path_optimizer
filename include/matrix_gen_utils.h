#ifndef MATRIX_GEN_UTILS_H_
#define MATRIX_GEN_UTILS_H_

/*
 * Type definitions
 */
typedef int c_int;
typedef double c_float;

enum MatFormat {RowMajor, ColMajor};

typedef struct {
  c_int    nzmax; ///< maximum number of entries.
  c_int    m;     ///< number of rows
  c_int    n;     ///< number of columns
  c_int   *p;     ///< column pointers (size n+1)
  c_int   *i;     ///< row indices, size nzmax starting from 0
  c_float *x;     ///< numerical values, size nzmax
  c_int    nz;    ///< # of entries in triplet matrix, -1 for csc
} csc;

typedef struct {
  c_int    nzmax; ///< maximum number of entries.
  c_int    m;     ///< number of rows
  c_int    n;     ///< number of columns
  c_int   *p;     ///< row pointers (size m+1)
  c_int   *i;     ///< column indices, size nzmax starting from 0
  c_float *x;     ///< numerical values, size nzmax
  c_int    nz;    ///< # of entries in triplet matrix, -1 for csc
} csr;

/*
 * Parameter
 */
static c_float RAND_INTERVAL = (c_float)10.0; // values of matrix in [-RAND_INTERVAL/2, +RAND_INTERVAL/2]
static c_float ROUND_TO_ZERO_BELOW = (c_float)1.0; // e.g. ROUND_TO_ZERO_BELOW=RAND_INTERVAL/4 gives 50% zero elements

/*
 * Dense Matrices
 */

// Allocates a matrix and fills it using (rand()/RAND_MAX - 0.5) * RAND_INTERVAL
c_float * random_dense_matrix(const c_int nrows, const c_int ncols, enum MatFormat format);

// Allocates a positive semi-definite matrix by squaring a matrix from random_dense_matrix()
c_float * random_dense_psd_matrix(const c_int nrows, enum MatFormat format);

// Sets all elements of a dense matrix with (element < ROUND_TO_ZERO_BELOW) to zero
void sparsify_dense_matrix(c_float * in_matrix, const c_int nrows, const c_int ncols);

// Print dense matrix to console
void print_dense_matrix(const c_float * in_matrix, const c_int nrows, const c_int ncols, const enum MatFormat format, const char *name);


/*
 * Compressed Sparse Row Format
 */

// Create a matrix in Compressed Sparse Row format using a dense matrix
csr * dense_to_csr_matrix(c_float * in_matrix, const c_int nrows, const c_int ncols, const enum MatFormat format);

// Create a dense matrix using a a CSR matrix
c_float * csr_to_dense_matrix(const csr * in_csr, const enum MatFormat out_format);

// Print CSR matrix to console
void print_csr_matrix(const csr * in_csr, const char *name);


/*
 * Compressed Sparse Column Format
 */

// Create a matrix in Compressed Sparse Column format using a dense matrix
csc * dense_to_csc_matrix(c_float * in_matrix, const c_int nrows, const c_int ncols, const enum MatFormat format);

// Create a dense matrix using a a CSC matrix
c_float * csc_to_dense_matrix(const csc * in_csc, const enum MatFormat out_format);

// Print CSC matrix to console
void print_csc_matrix(const csc * in_csc, const char *name);

#endif
