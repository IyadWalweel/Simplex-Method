import numpy as np
from itertools import chain
import copy
from tabulate import tabulate

def simplex_Method(c, table, x, obj, lengths, eps, phase1 = False, phase2 = False, basic_phs = None, Lev_I = [], Ent = [],
                                                                         Ent_I = [], Levx_I = [], XEnt = [], Big_M = False, M = 100):
    lc, lbl, lbg, lbe = lengths[0], lengths[1], lengths[2], lengths[3]
    C = []         # List of original objective coefficients of the Entering variables (Useful in the Two-Phase Only!)
    basic_phs2 = np.zeros(lbl+lbg+lbe)
    Entering = Ent.copy()    # List of x's Entering variables (To check if one variable was entering and then leaving)
    EI = Ent_I.copy() 
    LI = Levx_I.copy()
    I = Lev_I.copy()          # Indices of x's Entering Variables
    X = []    # List of Solutions 
    Y = []    # List of Dual Solutions
    F = []    # List of Function Values
    ite = 0
    # Xs = list(table[0][1:lc+1])
    Xs = XEnt.copy()
    # lc = len(c)
#     X = np.zeros((lc))
    # print("Entering = ", Entering)
    # print("EI = ", EI)
    # print("I = ", I)

    # print(tabulate(table))

#     O = np.zeros(lc+lbl+(2*lbg)+lbe + 1)
#     O[:lc] = c
    # O = list(table[1])
#     print("O = ",O)
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

    # else:
    #     if (obj == 'Max'):
    #         table[1][1:] *= -1          # Multiplying the Objective row by (-1) in case of Maximizing it 
    #         O = list(table[1])
    #         F.append(- O[-1])           # Initial Objective Function Value in case of Maximization  
    #     else: 
    #         O = list(table[1])
    #         print('O = ', O)
    #         F.append( O[-1])            # Initial Objective Function Value in case of Minimization
#     print(table[2:][-(lbl+lbg+lbe+1):-1])

#     x = np.zeros(lc)
#     z = x
    X.append(x.copy())              # Adding the Initial Solution (x) to the list of solutions (X)
#     print('X = ', X)
#     print(tabulate(table))
    OPI = table[2:, -(lbl+lbg+lbe+1):-1]      # Optimal Primal Inverse (OPI) - the initial matrix
#     print(tabulate(table))
    # print('OPI = ', OPI)
    if phase2:
        OOC = basic_phs
    else:
        OOC = O[-(lbl+lbg+lbe+1):-1]            # Original Objective Coefficients (OOC) of optimal primal basic variables - inital vector
#     print(tabulate(table))
#     print('type is', type(OOC))
    # print('OOC = ', OOC)
#     print("O = ", O)
#     print('OOC[0] = ', OOC[0])
#     ind1 = 0
#     ind2 = 3
#     OOC[ind1] = O[ind2]
#     print("OOC = ", OOC)
#     print("O = ", O)
#     print(OOC.shape)
    y = OOC@OPI                         # Dual Solution 
    Y.append(y)                         # Adding the inital dual solution (y) to the list of dual solutions (Y)
    # print("x = ", x)
    # print("y = ", y)
    
    ma = np.max(table[1][1:-1])              # Finding the maximum element in the objective row 
    # print('max = ', ma)
    ind = np.argmax(table[1][1:-1]) + 1      # Finding the index of this max. element 
    # print('max index = ', ind)
    while (ma > eps):                      # Stopping Criterion 
        entering = table[0][ind]             # Determining the Entering variable 
        Entering.append(entering)         # x's Entering variables

            # print(Entering)
        # print("EI = ", EI)
            # C[ite] = O[ind]                      # Original Objective coefficient of the entering variable (Two-Phase)
#         print(C)
#         print('entering is ', entering)
        L = []                                # List of Ratios to determing the leaving variable 
        for i in range(2,2+lbl+lbg+lbe):
            if (table[i][ind]>0):
                L.append(table[i][-1]/table[i][ind])
            else:
                L.append(float('inf'))
        # print('L = ', L)
        mi = min(L)
        mi_ind = L.index(mi) + 2             # Index of the leaving variable
        leaving = table[mi_ind][0]           # Leaving variable
        I.append(mi_ind)
        C.append(-O[ind])

        if (ind <= len(x)):
            EI.append(ind - 1)                # Indices of x's Entering variables
            Xs.append(entering)
            LI.append(mi_ind)
        
        basic_phs2[mi_ind-2] = O[ind]

        # print(leaving)
        # Leaving.append(leaving)
        # print('leav = ', leaving)
        # print('ent = ', entering)
        # I.append(mi_ind)                     # List of indices of the leaving variables 

        # print('leaving is ', leaving)
        # print('min index = ', mi_ind)
        table[mi_ind][0] = entering          
        # print(tabulate(table)) 
        pivot = table[mi_ind][ind]           # Pivot 
        # print('pivot is', pivot)
        table[mi_ind][1:] /= pivot           # Deviding by the Pivot 
        # print(tabulate(table))
        concatenated = chain(range(1,mi_ind), range(mi_ind+1, lbl+lbg+lbe+2))
        for j in concatenated:
            table[j][1:] += -table[j][ind]*table[mi_ind][1:]       
        # print(tabulate(table))
        if (obj == 'Min'):
            F.append(table[1][-1])
        else: 
            F.append(-table[1][-1])
        
        # print('entering is: ', entering)
        # print('leaving is: ', leaving)
        # print('Entering is: ', Entering)
        # print('EI = ', EI)
        # print('I = ', I)
        # print('Xs is: ', Xs)
        if leaving in Xs:
            k = Xs.index(leaving)
            x[EI[k]] = 0
            EI.remove(EI[k])
            LI.remove(LI[k])
            I.remove(I[k])
            if phase1:
                C.remove(C[k])
            Entering.remove(Entering[k])
            Xs.remove(Xs[k])
        elif leaving in Entering:
            k = Entering.index(leaving)
            # print('k = ', k)
            Entering.remove(Entering[k])
            # EI.remove(EI[k])
            I.remove(I[k])
            if phase1:
                C.remove(C[k])
            # print('EI New = ', EI)
            # print("I New = ", I)
            # if leaving in Xs:
            #     r = Xs.index(leaving)
            #     x[r] = 0                    # If leaving was entering make it's coordinate in the solution to be zero
            # I.append(mi_ind)
        #     C[k] = -O[ind]
        # else:
        #     # I.append(mi_ind)
        #     C[mi_ind-2] = -O[ind]
        # if entering in Xs:
        #     if leaving in Entering:
        #         k = Entering.index(leaving)
        #         C[k] = -O[ind]
        #     else:
        #         C[ite] = -O[ind]
        # print('I = ', I)
        # print('C = ', C)
        
        for i in range(len(EI)):
            x[EI[i]] = table[LI[i]][-1]         # Finding the solution 
        X.append(x.copy())
#         X = np.append(X,x)
#         if entering in table[0][1:lc+1]:
#             EI.append(ind - 1)
#             X[ind - 1] = table[mi_ind][-1]
#         print('x = ', x)
#         print('X = ', X)
        OPI = table[2:, -(lbl+lbg+lbe+1):-1]
        # print("OPI = ", OPI)
#         print(tabulate(table))
#         k = mi_ind - 2
#         print('k = ', k)
#         print('ind = ', ind)
        # print('O = ', O)
#         print('O[ind] = ', O[ind])
        OOC[mi_ind - 2] = O[ind]
        # print("OOC = ", OOC)
#         print(tabulate(table))
        y = OOC@OPI
        Y.append(y)
        # print("x = ", x)
        # print("y = ", y)
        # print('****************************')
        # print('ite = ', ite)
        ite += 1                                  # New Iteration
        ma = np.max(table[1][1:-1])
        # print('max = ', ma)
        ind = np.argmax(table[1][1:-1]) + 1
        # print('max index = ', ind)
        # print(tabulate(table))
    # print("iterations = ", ite)
    
    if phase1:
        print(table[1][-1])
        return table, X[-1], C, basic_phs2, I, Entering, EI, LI, Xs
    else:
        # print(len(EI))
        # print('EI = ', EI)
        # print(len(LI))
        # print("LI = ", LI)
        # print("X = ", X[-1])
        X = np.array(X)
        Y = np.array(Y)
        res = {
            'Optimal Value': F[-1],
            'Optimal Solution': X[-1],
            'Associated Dual Solution': Y[-1],
            'Number of Iterations': ite
        }
#         print(tabulate(res))
        return res, table, X, Y, F, ite
 