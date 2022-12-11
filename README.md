# LuDecomposition

Introduction
Often in computer science, it is important to solve a linear equation of the form A*x = B, with A being a matrix, x being a vector, and B being the solution vector. Rather than directly solving by matrix multiplication, A can be broken into two matrices L and U which is Lower Matrix and Upper Matrix. Lower Matrix which has non-zero values in a lower diagonal region of the matrix and the upper matrix which has non-zero values in the upper diagonal region of the matrix. 

The reason for doing this is to get the L matrix and U matrix. Then we can solve for y by using L*y = b (forward substitution), and then use U*x = y to solve for the vector x(background substitution). So with this method, we can solve for many vectors. It can also be used to find the determinant of a matrix. After finding L and U, it can now be used to get the product of their diagonal indexes, which can be used to find the determinant of A. 

The goal, however, is not just to write a regular, sequential algorithm that computes the L and U matrices for any non-invertible matrix A. That is just the first step. For the next three steps are to implement three parallel programs based on the sequential programs we create in step 1. Those three programs are Cuda (GPU-based), OpenMp (shared memory), MPI (distributed memory). For the original because LU decomposition has a format O(n^3), so for large matrix sizes, running them sequentially is not feasible. With the three parallel programs that can run parallel would make larger matrix sizes more practical to solve.
