#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix_gen_utils.h"

#define UNIT_TESTS

//#ifdef UNIT_TESTS
//extern void run_unit_tests();
//#endif

int test(void) {

	/* Example */
	const c_int nrows = 4;
	const c_int ncols = 5;

	// Create a random positive semi-definite matrix
	c_float * dense_matrix = random_dense_matrix(nrows, ncols, RowMajor);
	print_dense_matrix(dense_matrix, nrows, ncols, RowMajor, "dense matrix");

	// Sparsify matrix
	sparsify_dense_matrix(dense_matrix, nrows, ncols);
	print_dense_matrix(dense_matrix, nrows, ncols, RowMajor, "sparsified dense matrix");

	// Create a matrix in compressed sparse column format
	csc * csc_matrix = dense_to_csc_matrix(dense_matrix, nrows, ncols, RowMajor);
	print_csc_matrix(csc_matrix, "compressed sparse column matrix");

	// Transform back to dense matrix
	c_float * other_dense_matrix = csc_to_dense_matrix(csc_matrix, ColMajor);
	print_dense_matrix(other_dense_matrix, nrows, ncols, ColMajor, "other dense matrix");

	free(dense_matrix);
	free(other_dense_matrix);
	free(csc_matrix->x);
	free(csc_matrix->i);
	free(csc_matrix->p);
	free(csc_matrix);



	return 0;
}
