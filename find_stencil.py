
""" finds coefficient for symmetric stencil for M points of each side of midpoint """

from sympy.matrices import *
from sympy import factorial
from scipy.linalg import solve
import numpy as np

M = 3;
N = 2*M+1;

A = np.zeros((N,N))
d = np.zeros(N)

d[1] = 1;           # diff order 1

for i in range(N):
    for j in range(M+1):
         A[i,j] = A[i,-(j+1)] = (M-j)**i/(float)(factorial(i))   # values of matrix
         A[i,j] *= (-1)**(i*(j+1))                                     # alternating sign

a = np.linalg.solve(A, d)
print a
