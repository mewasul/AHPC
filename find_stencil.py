
""" finds coefficient for symmetric stencil for M points of each side of midpoint """

from sympy.matrices import *
from sympy import factorial

M = 3;
N = 2*M+1;

A = zeros(N,N)
d = zeros(N,1)

d[1] = 1;           # diff order 1

for i in range(N):
    for j in range(M+1):
         A[i,j] = A[i,-(j+1)] = (M-j)**i/(factorial(i))   # values of matrix
         A[i,j] *= (-1)**(i)                              # alternating sign

a = A.solve(d)

print a
