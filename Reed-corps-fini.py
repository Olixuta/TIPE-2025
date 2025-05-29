#Reed-Solomon

import numpy as np
from scipy.interpolate import lagrange
import struct
import numpy as np
from random import random

#inverse = pow(a, -1, 131)

def evaluation_li(i,x_values,eval,modulo=131):
    yi = 1 
    n=len(x_values)
    for j in range(n):
        if j!=i:
            yi = (yi * ((eval-x_values[j])%modulo) * (pow(x_values[i]-x_values[j],-1,modulo)%modulo))%modulo  # non inveerse modulo 7
    return yi % modulo
    
def lagrange(x,y,eval,modulo=131):
    assert len(x)==len(y)
    p_eval = 0 
    n = len(y)
    for i in range(n):
        p_eval =( p_eval + (y[i] % modulo)* evaluation_li(i,x,eval,modulo) ) % modulo 
    return p_eval
    
    
    

def poly_add(p1, p2, modulo):
    length = max(len(p1), len(p2))
    result = [0] * length
    for i in range(length):
        a = p1[i] if i < len(p1) else 0
        b = p2[i] if i < len(p2) else 0
        result[i] = (a + b) % modulo
    return result

def poly_mul(p1, p2, modulo):
    result = [0] * (len(p1) + len(p2) - 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            result[i + j] = (result[i + j] + p1[i] * p2[j]) % modulo
    return result

def lagrange_poly(x, y, modulo):
    n = len(x)
    final_poly = [0]

    for i in range(n):
        # Construct L_i(x)
        li_poly = [1]
        denom = 1
        for j in range(n):
            if i != j:
                # Multiply by (x - xj)
                li_poly = poly_mul(li_poly, [-x[j] % modulo, 1], modulo)
                # Multiply denominator (xi - xj)
                denom = (denom * (x[i] - x[j])) % modulo

        denom_inv = pow(denom, -1, modulo)
        # Multiply L_i(x) by yi * denom_inv
        li_scaled = [(coef * y[i] * denom_inv) % modulo for coef in li_poly]

        # Add to the final polynomial
        final_poly = poly_add(final_poly, li_scaled, modulo)

    return final_poly
