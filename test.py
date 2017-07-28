import numpy as np
import numpy.linalg as linalg
from sympy import *



def Hamilton(e,J):
    H=np.diag(e)
    H[0,len(e)-1]=J
    H[len(e)-1,0]=J
    a=J
    print('a',a)
    for index, value in enumerate(e):
            if(index < len(e)-1):
                H[index+1,index]=H[0,len(e)-1]
                H[index,index+1]=J
    return H
#
# def eigenwerte(M):
#     eigenValues,eigenVectors = linalg.eig(M)
#     idx = eigenValues.argsort()[::-1]
#     print(idx)
#     eigenValues = eigenValues[idx]
#     eigenVectors = eigenVectors[:,idx]
#     return eigenValues, eigenVectors

J , a= symbols('J , a')

Test_H=Hamilton([-a,a,-a,a],J)
Test_H_matrix=Matrix(Test_H)
print(Test_H)


#print('eigenwerte und vektoren nach numpy',eigenwerte(Test_H))


print('eigenwerte und eigenvektoren nach sympy',Test_H_matrix.eigenvects()[0])
print('eigenwerte und eigenvektoren nach sympy',Test_H_matrix.eigenvects()[1])
print('eigenwerte und eigenvektoren nach sympy',Test_H_matrix.eigenvects()[2])
print('eigenwerte und eigenvektoren nach sympy',Test_H_matrix.eigenvects()[3])
