import numpy as np
from .Simplex_Method import simplex_Method
from tabulate import tabulate

def phase1(c, table, x, obj, lengths):
    lc, lbl, lbg, lbe = lengths[0], lengths[1], lengths[2], lengths[3]
    for i in range(1,lbg+lbe+1):
        table[1][1:] += table[-i][1:] 
    # print(tabulate(table))
    table, x, C, basic_phs, I, Entering, EI, Xs = simplex_Method(c, table, x, obj, lengths, True)
#     print(tabulate(table))
    
    
    return table, x, C, basic_phs, I, Entering, EI, Xs

def phase2(c, table, x, obj, C, basic_phs, I, lengths, Entering, EI, Xs):
    print("Phase2:- ")
    lc, lbl, lbg, lbe = lengths[0], lengths[1], lengths[2], lengths[3]
    # print('C = ',C)
    table[1][-(lbg+lbe+1):-1] = -float('inf')
    # table = np.delete(table, np.s_[-(lbg+lbe+1):-1], 1)
    table[1][0] = 'z' 
    if (obj == 'Min'):       
        table[1][1:lc+1] = -c
    else:
        table[1][1:lc+1] = c
    table[1][-1] = 0
#     print(tabulate(table))
    for i in range(len(I)):
        table[1][1:] += C[i]*table[I[i]][1:]
#     print(tabulate(table))
    result, table, X, Y, F, ite = simplex_Method(c, table, x, obj, lengths, phase2 = True, basic_phs=basic_phs,
                                                                                 Lev_I = I, Ent = Entering, Ent_I = EI, XEnt= Xs)
#     print(tabulate(table))
    
    return result, table, X, Y, F, ite