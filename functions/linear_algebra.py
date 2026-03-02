import numpy as np

def gauss_jordan_elimination(matrix_A:np.ndarray,
                             vector_b:np.ndarray
                             )->tuple[np.ndarray,np.ndarray]:

    # b is the solution
    #could include the identity for the inverse
    A = matrix_A.copy()
    b = vector_b.copy()

    if len(A.T) != len(A):
        raise Warning('Matrix is not square')
    for i in range(len(A.T)):
        #i is column index
        # A[:,i] is column
        largest_index = None
        pivot = 0

        for j in range(i,len(A)):

            #j is from column index to row index max
            # loop over rows from row i
            # find max in column i
            # so A[j,i]
            if abs(A[j,i]) > abs(pivot):
                largest_index = j
                pivot = A[j,i]

        if abs(pivot) > 0:

            #swap then normalise
            A[[i, largest_index]]   = A[[largest_index, i]]
            A[i,:]                  = A[i,:]/pivot

            b[[i,largest_index]]    = b[[largest_index,i]]
            b[i]                    = b[i]/pivot

            #reduce above and below diagonal to 0
            # need to make a copy first for the b vector

            above_terms = A[:i, i].copy()
            below_terms = A[i+1:, i].copy()
            A[:i, :]    -= above_terms[:, None]     * A[i, :]
            A[i+1:, :]  -= below_terms[:, None]     * A[i, :]

            b[:i]       -= above_terms              * b[i]
            b[i+1:]     -= below_terms              * b[i]

        else:
            raise Warning('Matrix is not singular')

    return b

def crouts_algorithm(matrix_A):


    A = matrix_A.copy()
    B = np.zeros((len(A),len(A.T)))
    if len(A.T) != len(A):
        raise Warning('Matrix is not square')

    for i in range(len(A)):
        # for j in range(len(A.T)):
        #     if i == 0:
        #         B[i,j] = A[i,j]
        for j in range(i,A.shape[0]):
            B[j,i] = A[j,i]-np.sum(B[j,:i]*B[:i,i])
        for j in range(i+1,A.shape[0]):
            B[i,j] = (1/B[i,i])*(A[i,j]-np.sum(B[i,:i]*B[:i,j]))

    #B is the LU decomp
    return B

def crouts_algorithm_first(matrix_A):
    A = matrix_A.copy()
    B = np.zeros((len(A),len(A.T)))
    if len(A.T) != len(A):
        raise Warning('Matrix is not square')

    for i in range(len(A)):
        B[i, i] = 1
        for j in range(len(A.T)):
            if i == 0:
                B[i,j] = A[i,j]
            else:
                # checked the result of this against this link, so I think its right
                # www.emathhelp.net/calculators/linear-algebra/lu-decomposition-calculator
                if i <= j:
                    B[i,j] = A[i,j] - np.sum(B[i,:j]*B[:j,j])
                if i > j:
                    B[i,j] = (1/B[j,j])*(A[i,j]-np.sum(B[i,:j]*B[:j,j]))

    #B is the LU decomp
    return B

def forward_sub(matrix_A,y):
    z = np.zeros(matrix_A.shape[0])
    for i in range(matrix_A.shape[0]):
        z[i] = (1/matrix_A[i,i])*(y[i]-matrix_A[i,:i]@z[:i])
    return z

def backward_sub(matrix_A,z):
    c = np.zeros(matrix_A.shape[0])
    for i in range(matrix_A.shape[0]-1,-1,-1):
        c[i] = (z[i] - (matrix_A[i,i+1:]@c[i+1:]))
    return c

#todo
# def gaussian_with_back_sub():
# def LU_decomp():
# def SVD():

