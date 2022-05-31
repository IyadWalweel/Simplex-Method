import numpy as np

def titles(lc, lbl = 0, lbg = 0, lbe = 0):
    M = ["" for x in range(lc+lbl+(2*lbg)+lbe + 2)]
    M[0] = "Basic"
    M[-1] = "Solution"
    for i in range(1,lc+1):
        M[i] = "x{}".format(i)
    for i in range(lc+1, lc+lbg+1):
        M[i] = "p{}".format(i - lc)
    for i in range(lc+lbg+1, lc+lbg+lbl+1):
        M[i] = "s{}".format(i - lc - lbg)
    for i in range(lc+lbg+lbl+1, lc+lbl+(2*lbg)+lbe+1):
        M[i] = "R{}".format(i - (lc+lbg+lbl))
        
    R = ["" for x in range(lbl+lbg+lbe + 1)]
    R[0] = "z"
    for j in range(1, lbl+1):
        R[j] = "s{}".format(j)
    for j in range(lbl+1, lbg+lbe+lbl+1):
        R[j] = "R{}".format(j - lbl)
           
    return np.array([M]), np.array([R])
