#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "matrix_gen_utils.h"

c_float abs_val(const c_float in) {
	return (in >= 0.0) ? in : -in;
}

c_int arr_ind(const c_int i_col, const c_int i_row, const c_int nrows, const c_int ncols, const enum MatFormat format) {
	return (format == RowMajor) ? (i_col + i_row * ncols) : (i_row + i_col * nrows);
}

void dense_matrix_mpy(const c_float * in_mat1, const c_int nrows1, const c_int ncols1, const enum MatFormat format1,
					  const c_float * in_mat2, const c_int nrows2, const c_int ncols2, const enum MatFormat format2,
					  c_float * out_mat, const enum MatFormat format3) {
	c_int i_row, i_col, i_col2, ind1, ind2, ind3;
	c_float res;
	assert(ncols1 == nrows2);
	for (i_row = 0; i_row < nrows1; i_row++) {
		for (i_col2 = 0; i_col2 < ncols2; i_col2++) {
			res = 0.0;
			for (i_col = 0; i_col < ncols1; i_col++) {
				ind1 = arr_ind(i_col, i_row, nrows1, ncols1, format1);
				ind2 = arr_ind(i_col2, i_col, nrows2, ncols2, format2);
				res += in_mat1[ind1] * in_mat2[ind2];
			}
			ind3 = arr_ind(i_col2, i_row, nrows1, ncols2, format3);
			out_mat[ind3] = res;
		}
	}
}


/*
 * Copies a dense matrix to another dense matrix with the same MatFormat
 */
void copy_matrix(const c_float * in_matrix, c_float * out_matrix, const c_int nrows, const c_int ncols) {
	int i_row, i_col;
	for (i_row = 0; i_row < nrows; i_row++) {
		for (i_col = 0; i_col < ncols; i_col++) {
			out_matrix[i_col + i_row * ncols] = in_matrix[i_col + i_row * ncols];
		}
	}
}


/*
 * Transposes a dense matrix
 */
void transpose_matrix(c_float * in_out_matrix, const c_int nrows, const c_int ncols, const enum MatFormat format) {
	int i_row, i_col, ind_copy, ind;
	c_float * tmp_matrix = (c_float *)malloc(sizeof(c_float) * nrows * ncols);
	copy_matrix(in_out_matrix, tmp_matrix, nrows, ncols);
	for (i_row = 0; i_row < nrows; i_row++) {
		for (i_col = 0; i_col < ncols; i_col++) {
			ind_copy = arr_ind(i_col, i_row, nrows, ncols, format);
			ind = arr_ind(i_row, i_col, ncols, nrows, format);
			in_out_matrix[ind] = tmp_matrix[ind_copy];
		}
	}
	free(tmp_matrix);
}

/*
 * Allocates a matrix and fills it using (rand()/RAND_MAX - 0.5) * RAND_INTERVAL
 */
c_float * random_dense_matrix(const c_int nrows, const c_int ncols, enum MatFormat format) {
	int i_row, i_col;
	c_float rand_value;
	c_float * rand_matrix = (c_float *)malloc(nrows * ncols * sizeof(c_float));
	for (i_row = 0; i_row < nrows; i_row++) {
		for (i_col = 0; i_col < ncols; i_col++) {
			rand_value = (((c_float)rand() / ((c_float)RAND_MAX)) - 0.5) * RAND_INTERVAL;
			rand_matrix[arr_ind(i_col, i_row, nrows, ncols, format)] = rand_value;
		}
	}
	return rand_matrix;
}


/*
 * Allocates a matrix and fills it using (rand()/RAND_MAX - 0.5) * RAND_INTERVAL
 */
c_float * random_dense_psd_matrix(const c_int nrows, enum MatFormat format) {
	const int ncols = nrows;
	c_float * rand_matrix = random_dense_matrix(nrows, ncols, format);

	c_float * rand_matrix_copy = (c_float *)malloc(nrows * ncols * sizeof(c_float));
	c_float * rand_matrix_out = (c_float *)malloc(nrows * ncols * sizeof(c_float));

	copy_matrix(rand_matrix, rand_matrix_copy, nrows, ncols);
	transpose_matrix(rand_matrix_copy, nrows, ncols, format);

	dense_matrix_mpy(rand_matrix,      nrows, ncols, format,
					 rand_matrix_copy, ncols, nrows, format,
					 rand_matrix_out, format);

	free(rand_matrix);
	free(rand_matrix_copy);
	return rand_matrix_out;
}


/*
 * Sets all elements < ROUND_TO_ZERO_BELOW of a matrix to zero.
 */
void sparsify_dense_matrix(c_float * in_matrix, const c_int nrows, const c_int ncols) {
	int i_row, i_col;
	for (i_row = 0; i_row < nrows; i_row++) {
		for (i_col = 0; i_col < ncols; i_col++) {
			if (abs_val(in_matrix[i_col + i_row * ncols]) < ROUND_TO_ZERO_BELOW) {
				in_matrix[i_col + i_row * ncols] = 0.0;
			}
		}
	}
}


/*
 * Get number of nonzero elements in a dense matrix
 */
c_int nonzero_elements(const c_float * in_matrix, const c_int nrows, const c_int ncols) {
	int i_row, i_col;
	c_int num_nonzero = 0;
	for (i_row = 0; i_row < nrows; i_row++) {
		for (i_col = 0; i_col < ncols; i_col++) {
			if (in_matrix[i_col + i_row * ncols] != 0.0) {
				num_nonzero++;
			}
		}
	}
	return num_nonzero;
}


/*
 * Copy and change matrix format
 */
void change_major_format(const c_float * in_matrix, const enum MatFormat in_format,
					 	 c_float * out_matrix, const enum MatFormat out_format,
					 	 const c_int nrows, const c_int ncols) {
	int i_row, i_col;
	for (i_col = 0; i_col < ncols; i_col++) {
		for (i_row = 0; i_row < nrows; i_row++) {
			out_matrix[arr_ind(i_col, i_row, nrows, ncols, out_format)] =
					in_matrix[arr_ind(i_col, i_row, nrows, ncols, in_format)];
		}
	}
}


/*
 * Dense to CSR matrix
 */
csr * dense_to_csr_matrix(c_float * in_matrix, const c_int nrows, const c_int ncols, const enum MatFormat format) {
	c_int i_row, i_col, ind_mat, ind_val = 0, num_row_nnz = 0;
	c_float *values;
	c_int *columns, *row_nnz;
	const c_int nnz = nonzero_elements(in_matrix, nrows, ncols);
	values = (c_float *)malloc(sizeof(c_float) * nnz);
	columns = (c_int *)malloc(sizeof(c_int) * nnz);
	row_nnz = (c_int *)malloc(sizeof(c_int) * (nrows + 1));

	// Fill values
	row_nnz[0] = (c_int)0;
	for (i_row = 0; i_row < nrows; i_row++) {
		num_row_nnz = 0;
		for (i_col = 0; i_col < ncols; i_col++) {
			ind_mat = arr_ind(i_col, i_row, nrows, ncols, format);
			if (in_matrix[ind_mat] != 0.0) {
				values[ind_val] = in_matrix[ind_mat];
				columns[ind_val] = i_col;
				ind_val++;
				num_row_nnz++;
			}
		}
		row_nnz[i_row + 1] = row_nnz[i_row] + num_row_nnz;
	}

	// Create CSR structure
	csr * csr_matrix = (csr *)malloc(sizeof(csr));
	csr_matrix->nzmax = nnz;
	csr_matrix->nz = nnz;
	csr_matrix->m = nrows;
	csr_matrix->n = ncols;
	csr_matrix->p = row_nnz;
	csr_matrix->x = values;
	csr_matrix->i = columns;

	return csr_matrix;
}

/*
 * CSR to dense matrix
 */
c_float * csr_to_dense_matrix(const csr * in_csr, const enum MatFormat out_format) {
	c_int i_row, i_col, nnz_in_row, i_val = 0, i_nnz;
	const c_int nrows = in_csr->m;
	const c_int ncols = in_csr->n;
	c_float * out_matrix = (c_float *)calloc(nrows * ncols, sizeof(c_float));
	c_int * row_nnz = in_csr->p;
	c_int * columns = in_csr->i;
	c_float * values = in_csr->x;

	for (i_row = 0; i_row < nrows; i_row++) {
		nnz_in_row = row_nnz[i_row + 1] - row_nnz[i_row];
		if (nnz_in_row > 0) {
			for (i_nnz = 0; i_nnz < nnz_in_row; i_nnz++) {
				i_col = columns[i_val];
				out_matrix[arr_ind(i_col, i_row, nrows, ncols, out_format)] = values[i_val];
				i_val++;
			}
		}
	}

	return out_matrix;
}


/*
 * Dense to CSC matrix
 */
csc * dense_to_csc_matrix(c_float * in_matrix, const c_int nrows, const c_int ncols, const enum MatFormat format) {
	c_int i_row, i_col, ind_mat, ind_val = 0, num_col_nnz = 0;
	c_float *values;
	c_int *rows, *col_nnz;
	const c_int nnz = nonzero_elements(in_matrix, nrows, ncols);
	values = (c_float *)malloc(sizeof(c_float) * nnz);
	rows = (c_int *)malloc(sizeof(c_int) * nnz);
	col_nnz = (c_int *)malloc(sizeof(c_int) * (ncols + 1));

	// Fill values
	col_nnz[0] = (c_int)0;
	for (i_col = 0; i_col < ncols; i_col++) {
		num_col_nnz = 0;
		for (i_row = 0; i_row < nrows; i_row++) {
			ind_mat = arr_ind(i_col, i_row, nrows, ncols, format);
			if (in_matrix[ind_mat] != 0.0) {
				values[ind_val] = in_matrix[ind_mat];
				rows[ind_val] = i_row;
				ind_val++;
				num_col_nnz++;
			}
		}
		col_nnz[i_col + 1] = col_nnz[i_col] + num_col_nnz;
	}

	// Create CSR structure
	csc * csc_matrix = (csc *)malloc(sizeof(csc));
	csc_matrix->nzmax = nnz;
	csc_matrix->nz = nnz;
	csc_matrix->m = nrows;
	csc_matrix->n = ncols;
	csc_matrix->p = col_nnz;
	csc_matrix->x = values;
	csc_matrix->i = rows;

	return csc_matrix;
}

/*
 * CSC to dense matrix
 */
c_float * csc_to_dense_matrix(const csc * in_csc, const enum MatFormat out_format) {
	c_int i_row, i_col, nnz_in_col, i_val = 0, i_nnz;
	const c_int nrows = in_csc->m;
	const c_int ncols = in_csc->n;
	c_float * out_matrix = (c_float *)calloc(nrows * ncols, sizeof(c_float));
	c_int * col_nnz = in_csc->p;
	c_int * rows = in_csc->i;
	c_float * values = in_csc->x;

	for (i_col = 0; i_col < ncols; i_col++) {
		nnz_in_col = col_nnz[i_col + 1] - col_nnz[i_col];
		if (col_nnz > 0) {
			for (i_nnz = 0; i_nnz < nnz_in_col; i_nnz++) {
				i_row = rows[i_val];
				out_matrix[arr_ind(i_col, i_row, nrows, ncols, out_format)] = values[i_val];
				i_val++;
			}
		}
	}

	return out_matrix;
}


void print_dense_matrix(const c_float * in_matrix, const c_int nrows, const c_int ncols, const enum MatFormat format, const char *name) {
	int i_col, i_row, ind;
	printf("%s =\n[", name);
	for (i_row = 0; i_row < nrows; i_row++) {
		if (i_row > 0) {
			printf(" ");
		}
		for (i_col = 0; i_col < ncols; i_col++) {
			ind = arr_ind(i_col, i_row, nrows, ncols, format);
			if ((i_col == ncols-1) && (i_row == nrows-1)) {
				if (in_matrix[ind] >= 0) {
					printf(" %.4f", in_matrix[ind]);
				} else {
					printf("%.4f", in_matrix[ind]);
				}
			} else {
				if (in_matrix[ind] >= 0) {
					printf(" %.4f,", in_matrix[ind]);
				} else {
					printf("%.4f,", in_matrix[ind]);
				}
			}
		}
		if (i_row == nrows - 1) {
			printf("]\n");
		} else {
			printf("\n");
		}
	}
}


void print_csc_matrix(const csc * in_csc, const char *name) {
	c_int i;
	c_int * col_nnz = in_csc->p;
	c_int * rows = in_csc->i;
	c_float * values = in_csc->x;
	printf("%s =\n{", name);
	printf("   values = [");
	for (i = 0; i < in_csc->nz; i++) {
		if (i == in_csc->nz - 1) {
			printf("%.2f]\n", values[i]);
		} else {
			printf("%.2f, ", values[i]);
		}
	}
	printf("      rows = [");
	for (i = 0; i < in_csc->nz; i++) {
		if (i == in_csc->nz - 1) {
			printf("%d]\n", rows[i]);
		} else {
			printf("%d, ", rows[i]);
		}
	}
	printf("   col_nnz = [");
	for (i = 0; i < in_csc->n + 1; i++) {
		if (i == in_csc->n) {
			printf("%d]\n", col_nnz[i]);
		} else {
			printf("%d, ", col_nnz[i]);
		}
	}
	printf("}\n");
}


void print_csr_matrix(const csr * in_csr, const char *name) {
	c_int i;
	c_int * row_nnz = in_csr->p;
	c_int * cols = in_csr->i;
	c_float * values = in_csr->x;
	printf("%s =\n{", name);
	printf("   values = [");
	for (i = 0; i < in_csr->nz; i++) {
		if (i == in_csr->nz - 1) {
			printf("%.2f]\n", values[i]);
		} else {
			printf("%.2f, ", values[i]);
		}
	}
	printf("      cols = [");
	for (i = 0; i < in_csr->nz; i++) {
		if (i == in_csr->nz - 1) {
			printf("%d]\n", cols[i]);
		} else {
			printf("%d, ", cols[i]);
		}
	}
	printf("   col_nnz = [");
	for (i = 0; i < in_csr->m + 1; i++) {
		if (i == in_csr->m) {
			printf("%d]\n", row_nnz[i]);
		} else {
			printf("%d, ", row_nnz[i]);
		}
	}
	printf("}\n");
}

