from Titles import titles
import numpy as np
from tabulate import tabulate

def std_Form(obj, c, Al = None, bl = None, Ag = None, bg = None, Ae = None, be = None, Show = False, Big_M = False, M = 100):
    lc = len(c)
    # M = 0                         # In case Big_M = False
    case = None                   # Determining the case of the problem: simple simplex, Big-M, Two Phase
    lengths = np.zeros(4, int)
    lengths[0] = lc
    if ((bg is None) and (be is None)):     # The Simplist case 
        case = 'slack'
        lbl = len(bl)
        Id = np.identity(lbl)              # Identity for the slack variables
        if (lbl == 1):                     # Make Al and bl two-dim. arrays 
            Al = np.array([Al])            
            bl = np.array([bl])
            T = np.concatenate((Al,Id,bl.T), axis = 1, dtype = 'O')
        else:
            bl = np.array([bl])
            T = np.concatenate([Al,Id,bl.T], axis = 1, dtype = 'O')
        O = np.zeros(lc+lbl+1) 
        if (obj == "Min"):
            O[:lc] = -c                   # Objective function row 
        else:
            O[:lc] = c
        T = np.concatenate(([O],T))
        ct, rt = titles(lc,lbl)        # Column's and row's titles
        T = np.concatenate((rt.T, T), axis = 1)
        T = np.concatenate((ct, T))
        if Show:
            print('Problem in Standard Form:')
            print(tabulate(T))
        lengths[1] = lbl
        return T, case, lengths
    else:
        B = [bl,bg,be]
        A = [Al,Ag,Ae]
        for i in range(len(B)):
            if (B[i] is not None):
                lengths[i+1] = len(B[i])
                if (lengths[i+1] == 1):
                    A[i] = np.array([A[i]])      # Make it two-dim. array
        A = list(filter(lambda ele:ele is not None, A))        # Deleting the None arrays
        B = list(filter(lambda ele:ele is not None, B))        # Deleting the None arrays
        lbl, lbg, lbe = lengths[1], lengths[2], lengths[3]
        A = np.concatenate(A, dtype = 'O')       # Concatenate Al, Ag, and Ae
        B = np.array([np.concatenate(B)])        # Concatenate bl, bg, and be
        Id = np.identity(lbl+lbg+lbe)            # Identity of all the slack, and artificial variables 
        if (bg is not None):                     # Adding the surplus variables
            P = np.zeros((lbl+lbg+lbe, lbg))
            for i in range(lbg):
                P[lbl+i][i] = -1
            T = np.concatenate((A,P,Id, B.T), axis = 1)
        else:
            T = np.concatenate((A, Id, B.T), axis = 1)
        ct, rt = titles(lc,lbl,lbg,lbe)                  # Column's and row's titles
        if (Big_M):
            # print("Std-Form -- Big M")
            O = np.zeros(lc+lbl+(2*lbg)+lbe + 1)   
            O[-(lbg+lbe+1):-1] = -M*np.ones(lbg+lbe)         # Objective function row with Big-M 
            # print(lbe)
            # print(lbg)
            # print('O in Big M is',O)  
            if (obj == "Min"):
                O[:lc] = -c
            else:
                O[:lc] = c
            # print("O",O)
            O = np.array([O])
        else: # Two-Phases Method
            O = np.zeros(lc+lbl+(2*lbg)+lbe + 1)
            O[-(lbg+lbe+1):-1] = -np.ones(lbg+lbe)           # Objective function row with Two-Phase
            O = np.array([O])
            rt[0][0] = 'r'
        T = np.concatenate((O,T))
        T = np.concatenate((rt.T, T),axis = 1)
        T = np.concatenate((ct, T))
        if Show:
            print('Problem in Standard Form:')
            print(tabulate(T))
        return T, case, lengths