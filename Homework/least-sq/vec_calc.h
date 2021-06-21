double gsl_vector_length(gsl_vector* vec);

void matrix_print(FILE* stream, const gsl_matrix *X);

void vector_print(FILE* stream, gsl_vector* vec);

void GS_decomp(gsl_matrix *A, gsl_matrix *B);

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x);

void QR_inverse(gsl_matrix *Q, gsl_matrix *B);
