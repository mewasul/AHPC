
""" finds coefficient for symmetric stencil for M points of each side of midpoint """

from sympy.matrices import *
from sympy import factorial
from scipy.linalg import solve

M = 3;
N = 2*M+1;

A = zeros(N,N)
d = zeros(N)

d[1] = 1;           # diff order 1

for i in range(N):
    for j in range(M+1):
         A[i,j] = A[i,-(j+1)] = (M-j)**i/(factorial(i))   # values of matrix
         A[i,j] *= (-1)**(i*(j+1))                                     # alternating sign


print A
"""
print det(A[:1,:1])
print det(A[:2,:2])
print det(A[:3,:3])
print det(A[:4,:4])
print det(A[:5,:5])
print det(A[:6,:6])
"""
