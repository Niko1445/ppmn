void CB_decomb(gsl_matrix* A, gsl_matrix *L);

void CB_Asolve(gsl_matrix *A, gsl_vector *b, gsl_vector *x);

void CB_solve(gsl_matrix *L, gsl_vector *b, gsl_vector *x);

void CB_inverse(gsl_matrix *A, gsl_matrix *B);

double CB_det(gsl_matrix *A);
