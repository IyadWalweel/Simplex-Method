import numpy as np
from Standard_Form import std_Form
from Simplex_Method import simplex_Method
from Big_M import Big_M_Method
from Two_Phase import *
from tabulate import tabulate

def simplex(c, Al = None, bl = None, Ag = None, bg = None, Ae = None, be = None, obj = 'Min', Show = False, Big_M = False, M=100):
    B = [bl,bg,be]
    S = ['bl','bg','be']
    for i in range(len(B)):
        if B[i] is not None:
            assert B[i].any()>=0, f'All the {S[i]} = {B[i]} entries should be non-negative!'
    table, case, lengths = std_Form(obj, c, Al, bl, Ag, bg, Ae, be, Show, Big_M, M)
    lc, lbl, lbg, lbe = lengths[0], lengths[1], lengths[2], lengths[3]
    if (case is not None):
        x = np.zeros(lc)
        result, table, X, Y, F, ite = simplex_Method(c, table, x, obj, lengths)
    elif (Big_M):
        x = np.zeros(lc)
        result, table, X, Y, F, ite = Big_M_Method(c, table, M, x, obj, lengths)
    else:
        x = np.zeros(lc)
        table, x_new, C, I = phase1(c, table, x, obj, lengths)
        result, table, X, Y, F, ite = phase2(c, table, x_new, obj, C, I, lengths)
    if Show:
        print("The Optimal Tableau: ")
        print(tabulate(table))
        
    return result, table, X, Y, F, ite