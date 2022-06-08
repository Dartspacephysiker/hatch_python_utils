import numpy as np

def can_geninverse_left(A):
    rankA = np.linalg.matrix_rank(A)
    return rankA == A.shape[1]

def can_geninverse_right(A):
    rankA = np.linalg.matrix_rank(A)
    return rankA == A.shape[0]

def geninverse_left(A):

    if can_geninverse_left(A):
        return np.linalg.inv(A.T@A)@(A.T)
    else:
        print("Generalized left-inverse of A not possible!")
        return A*0
    
def geninverse_right(A):

    if can_geninverse_right(A):
        return (A.T)@np.linalg.inv(A@A.T)
    else:
        print("Generalized right-inverse of A not possible!")
        return A*0

if __name__ == '__main__':

    A = np.array([[1,2,3,4],
                  [4,5,6,-7]])
    # A = np.array([[1,0,0,0],
    #               [0,1,0,0]])
    print(A)
    print("A   shape:",A.shape)
    print("A    rank:",np.linalg.matrix_rank(A))
    print("A^T  rank:",np.linalg.matrix_rank(A.T))
    
    
    AL = geninverse_left(A)
    AR = geninverse_right(A)
    print("")
    print("Gen. inverse left")
    print(AL)
    print("")
    print("Gen. inverse right")
    print(AR)
    
    
