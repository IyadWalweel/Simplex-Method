from .Simplex_Method import simplex_Method

def Big_M_Method(c, table, M, x, obj, lengths, eps):
    lc, lbl, lbg, lbe = lengths[0], lengths[1], lengths[2], lengths[3]
    for i in range(1,lbg+lbe+1):
        # if (obj == "Min"):
        table[1][1:] += M*table[-i][1:]
        # else:
        #     table[1][1:] -= M*table[-i][1:]
#     print(tabulate(table))
    result, table, X, Y, F, ite = simplex_Method(c, table, x, obj, lengths, eps, Big_M=True, M = M)
    
    return result, table, X, Y, F, ite