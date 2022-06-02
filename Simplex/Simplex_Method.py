import numpy as np
from itertools import chain
import copy
from tabulate import tabulate

def simplex_Method(c, table, x, obj, lengths, eps, phase1 = False, phase2 = False, basic_phs = None, Lev_I = [], Ent = [],
                                                                         Ent_I = [], Levx_I = [], XEnt = [], Big_M = False, M = 100):
    lc, lbl, lbg, lbe = lengths[0], lengths[1], lengths[2], lengths[3]
    basic_phs2 = np.zeros(lbl+lbg+lbe)
    Entering = Ent.copy()    # List of x's Entering variables (To check if one variable was entering and then leaving)
    EI = Ent_I.copy() 
    LI = Levx_I.copy()
    I = Lev_I.copy()          # Indices of x's Entering Variables
    X = []    # List of Solutions 
    Y = []    # List of Dual Solutions
    F = []    # List of Function Values
    ite = 0
    Xs = XEnt.copy()


    O = np.zeros(lc+lbl+(2*lbg)+lbe+2)      
    if (obj == 'Max'):
        F.append(- table[1][-1])
        O[1:lc+1] = c
        if Big_M:
            O[-(lbg+lbe+1):-1] = -M*np.ones(lbg+lbe)
    else:
        F.append( table[1][-1])
        O[1:lc+1] = -c
        if Big_M:
            O[-(lbg+lbe+1):-1] = -M*np.ones(lbg+lbe)


    X.append(x.copy())              # Adding the Initial Solution (x) to the list of solutions (X)
    OPI = table[2:, -(lbl+lbg+lbe+1):-1]      # Optimal Primal Inverse (OPI) - the initial matrix
    if phase2:
        OOC = basic_phs
    else:
        OOC = list(O[-(lbl+lbg+lbe+1):-1])            # Original Objective Coefficients (OOC) of optimal primal basic variables - inital vector

    y = OOC@OPI                         # Dual Solution 
    Y.append(y)                         # Adding the inital dual solution (y) to the list of dual solutions (Y)
    
    ma = np.max(table[1][1:-1])              # Finding the maximum element in the objective row 
    ind = np.argmax(table[1][1:-1]) + 1      # Finding the index of this max. element 
    while (ma > eps):                      # Stopping Criterion 
        entering = table[0][ind]             # Determining the Entering variable 
        Entering.append(entering)         # x's Entering variables

        L = []                                # List of Ratios to determing the leaving variable 
        for i in range(2,2+lbl+lbg+lbe):
            if (table[i][ind]>0):
                L.append(table[i][-1]/table[i][ind])
            else:
                L.append(float('inf'))
        mi = min(L)
        mi_ind = L.index(mi) + 2             # Index of the leaving variable
        leaving = table[mi_ind][0]           # Leaving variable
        I.append(mi_ind)

        if (ind <= len(x)):
            EI.append(ind - 1)                # Indices of x's Entering variables
            Xs.append(entering)
            LI.append(mi_ind)
        
        basic_phs2[mi_ind-2] = O[ind].copy()
        table[mi_ind][0] = entering          
        pivot = table[mi_ind][ind]           # Pivot 
        table[mi_ind][1:] /= pivot           # Deviding by the Pivot 
        concatenated = chain(range(1,mi_ind), range(mi_ind+1, lbl+lbg+lbe+2))
        for j in concatenated:
            table[j][1:] += -table[j][ind]*table[mi_ind][1:]       
        if (obj == 'Min'):
            F.append(table[1][-1])
        else: 
            F.append(-table[1][-1])

        if leaving in Xs:
            k1 = Xs.index(leaving)
            k2 = Entering.index(leaving)
            x[EI[k1]] = 0
            EI.remove(EI[k1])
            LI.remove(LI[k1])
            I.remove(I[k2])
            Entering.remove(Entering[k2])
            Xs.remove(Xs[k1])
        elif leaving in Entering:
            k = Entering.index(leaving)
            Entering.remove(Entering[k])
            I.remove(I[k])

        for i in range(len(EI)):
            x[EI[i]] = table[LI[i]][-1]         # Finding the solution 
        X.append(x.copy())

        OPI = table[2:, -(lbl+lbg+lbe+1):-1]
        OOC[mi_ind - 2] = O[ind].copy()
        y = OOC@OPI
        Y.append(y)
        ite += 1                                  # New Iteration
        ma = np.max(table[1][1:-1])
        ind = np.argmax(table[1][1:-1]) + 1
    
    if phase1:
        return table, X[-1], basic_phs2, I, Entering, EI, LI, Xs
    else:
        X = np.array(X)
        Y = np.array(Y)
        res = {
            'Optimal Value': F[-1],
            'Optimal Solution': X[-1],
            'Associated Dual Solution': Y[-1],
            'Number of Iterations': ite
        }
        return res, table, X, Y, F, ite
 