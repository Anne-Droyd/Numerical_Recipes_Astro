import numpy as np
from scipy.special.cython_special import log_wright_bessel


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


def crout(A):
    # Creating two L and U matrices filled with 0s and the same size as A

    L = np.zeros([A.shape[0],A.shape[0]])
    U = np.zeros([A.shape[0],A.shape[0]])
    n = len(A)

    # for-loop in order to set the j,j-th entry of U to 1

    for z in range(n):

        U[z, z] = 1

        # for-loop starting at L[j][j] in order to solve the j-th column of L

        for j in range(z, n):

            # Declaring a temporary L to store values and insert them later in the L matrix

            temporary_L = float(A[j, z])

            for k in range(z):
                temporary_L -= L[j, k] * U[k, z]

            L[j, z] = temporary_L

        # for-loop starting at U[j][j+1] in order to solve the j-th row of U

        for j in range(z + 1, n):

            # Declaring a temporary U to store values and insert them later in the U matrix

            temporary_U = float(A[z, j])

            for k in range(z):
                temporary_U -= L[z, k] * U[k, j]

            U[z, j] = temporary_U / L[z, z]

    # Returning the matrices L and U i.e. the A matrix decomposed using Crout's algorithm

    return (L, U)

def forward_substitution(L, b):
    y = np.full_like(b, 0)  # Creating the y vector the same size as the b vector

    for k in range(len(b)):

        y[k] = b[k]

        for i in range(k):
            y[k] = y[k] - (L[k, i] * y[i])

        y[k] = y[k] / L[k, k]  # Using forward substitution to find the value of y

    # Returning the y vector

    return y


def backward_substitution(U, y):
    x = np.full_like(y, 0)  # Creating the x vector the same size as the y vector

    for k in range(len(x), 0, -1):
        x[k - 1] = (y[k - 1] - np.dot(U[k - 1, k:], x[k:])) / U[
            k - 1, k - 1]  # Using backward substitution to find the value of x

    # Returning the solution vector x

    return x

def forward_sub(lower,y):
    #I think this is correct...
    #create a vector for transitioning between forward and backward substitution to find C
    #the choice of shape might matter here if the matrix is not square but ours is...
    z=np.zeros(lower.shape[0])
    for i in range(len(z)):
        #double checked but this seems correct...
        inverse=1/(lower[i][i])
        z[i]=(y[i]-(lower[i,:i]@z[:i]))*inverse

    return z

def backward_sub(upper,z):
    #double checked and this seems correct
    #create the c vector but not sure if this is to be vertical/horizontal
    c=np.zeros(upper.shape[0])

    #opposite of the last one have to go from bottom to top right to left
    #so need to start at the last item [len(c)-1] ie 19, to the first item ie 0, in negative steps
    #works better if we loop to the -1 item which doesn't make sense.... somehow goes to 0... not sure why
    #I think this is where it breaks and I got lucky in thinking it was correct

    for i in range(len(c)-1,-1,-1):
        #this is just 1 so no need
        inverse=1/upper[i][i]
        #double checked and this seems fine... need the +1 term
        c[i]=(z[i]-(upper[i,i+1:]@c[i+1:]))
    return c

# def forward_sub(matrix_A,y):
#     z = np.zeros(matrix_A.shape[0])
#     for i in range(matrix_A.shape[0]):
#         z[i] = (1/1)*(y[i]-matrix_A[i,:i]@z[:i])
#     return z
#
# def backward_sub(matrix_A,z):
#     c = np.zeros(matrix_A.shape[0])
#     for i in range(matrix_A.shape[0]-1,-1,-1):
#         c[i] = (1/matrix_A[i,i])*(z[i] - (matrix_A[i,i+1:]@c[i+1:]))
#     return c

#todo
# def gaussian_with_back_sub():
# def LU_decomp():
# def SVD():

