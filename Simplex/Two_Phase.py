import numpy as np
from .Simplex_Method import simplex_Method
from tabulate import tabulate

def phase1(c, table, x, obj, lengths, eps):
    lc, lbl, lbg, lbe = lengths[0], lengths[1], lengths[2], lengths[3]
    for i in range(1,lbg+lbe+1):
        table[1][1:] += table[-i][1:] 
    table, x, basic_phs, I, Entering, EI, LI, Xs = simplex_Method(c, table, x, obj, lengths, eps, True)    
    
    return table, x, basic_phs, I, Entering, EI, LI, Xs

##########

def phase2(c, table, x, obj, basic_phs, I, lengths, eps, Entering, EI, LI, Xs):
    print("Phase2:- ")
    lc, lbl, lbg, lbe = lengths[0], lengths[1], lengths[2], lengths[3]
    table[1][-(lbg+lbe+1):-1] = -float('inf')
    table[1][0] = 'z' 
    if (obj == 'Min'):       
        table[1][1:lc+1] = -c
    else:
        table[1][1:lc+1] = c
    table[1][-1] = 0
    for i in range(len(basic_phs)):
        table[1][1:] -= basic_phs[i]*table[i+2][1:]
    result, table, X, Y, F, ite = simplex_Method(c, table, x, obj, lengths, eps, phase2 = True, basic_phs=basic_phs, 
                                                        Lev_I = I, Ent = Entering, Ent_I = EI, Levx_I= LI, XEnt= Xs)    
    
    return result, table, X, Y, F, ite