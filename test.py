from Simplex_Solver import simplex as splx 
import numpy as np

# c = np.array([-7, -9, -18, -17])
# A = np.array([[2,4,6,7],[1,1,2,2],[1,2,3,3]])
# b = np.array([41,17,24])

c = np.array([5,12,4])
Al = np.array([1,2,1])
bl = np.array([10])
Ae = np.array([2,-1,3])
be = np.array([8])

result, table, X, Y, F, ite = splx(c,Al,bl,Ae=Ae,be=be,Show = True, Big_M=True, obj='Max')
print(result)

