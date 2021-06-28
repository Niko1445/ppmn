This is a implementation of the Cholesky-Banachiewicz (CB) decomposition algorithm using gsl vectors, matrices and BLAS support.
The matrix A is made to be an (n, n) upscaleable real symmetric positive definite matrix with a determinant not equal to 1,
by changing the parameter "n" from the top of the main a new A with size (n, n) can be made.
Other than the CB algorithm, a linear equation solver is implemented as well as a function for calculating the inverse matrix
and a function for calculating the determinant of a matrix, all using the CB decomposition algorithm, meaning they are limited
to real symmetric positive definite matrices.

Two solver functions are implemented, one takes matrix A, does the decomposition and solves, the other takes L and solves without a decomposition.
The latter is made such that if a function (such as CB_inverse) makes multiple calls to the solve function, the CB_decomp
won't be called for every solve but simply once. The former is simply for convenience.

checks.txt shows that the values calculated are correct and that unsurprisingly, gsl function are way faster than mine.

Optional self-evaluation out of 10: 8
