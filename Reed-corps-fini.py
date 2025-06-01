#Reed-Solomon

import numpy as np
from scipy.interpolate import lagrange
import struct
import numpy as np
from random import random


def repr_bin_str(c:str):
    """converti un charactère en sa représentation binaire

    Args:
        c (str): un str de taille 1 

    Returns:
        str: représentation binaire de c 
    """
    return format((ord(c)),'07b')


def codage_binaire(chaine:str):
    """converti une chaine en une représentation bianire 

    Args:
        chaine (str): 

    Returns:
        str: concaténation de la représentation binaire de chacun des caractère de la chaine en entrée 
    """
    c=''
    for elt in chaine:
        c=c + repr_bin_str(elt)
    return c




def decoupage_taille_p(message:str,p:int):# yes
    n = len(message)
    q= n//p
    r = n%p
    l = ['' for i in range(p+1)]
    for i in range(p+1):
        c=''
        for j in range(q):
            if i*q + j < n : 
                c=c+ message[i*q + j]
            #print(i*q+j)
        l[i]=c
    if r != 0 :
        #print((p)*q + q)
        l.append(message[(p)*q+q:])# le dernier est 
    return l 

m= 'This is a simple ASCII text '
n = len(m)
print(decoupage_taille_p(m,n))

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


def format_polynomial(coeffs, modulo):
    terms = []
    for i, coeff in enumerate(coeffs):
        coeff = coeff % modulo
        if coeff == 0:
            continue
        if i == 0:
            term = f"{coeff}"
        elif i == 1:
            term = f"{coeff}x" if coeff != 1 else "x"
        else:
            term = f"{coeff}x^{i}" if coeff != 1 else f"x^{i}"
        terms.append(term)
    
    if not terms:
        return "0"
    return " + ".join(terms[::-1]) + f"  (mod {modulo})"


x= [1,2,3,4]
y =[3,1,5,0]
print(lagrange(x,y,5,7))
p = lagrange_poly(x,y,7)
print(format_polynomial(p,7))

x= [2,3,4,5]
y =[1,5,0,6]
print(lagrange(x,y,1,7))
p=lagrange_poly(x,y,7)
print(format_polynomial(p,7))


# a modif 
def encodage_reed_solomon(message:str,redondance:int,p:int):
    n = len(message)
    l = decoupage_taille_p(message,n)
    ll=[ 0 for i in range(m+redondance)]
    m = len(m)
    x_points=[i for i in range(m)]
    y_points=[ord(message[i]) for i in range(n)]

    P=lagrange_poly(x_points,y_points,131)
    for r in range(redondance):
        #print(P(m+r))
        #assert(abs(P(m+r))<=2**15)
        l.append(P(m+r))

    return l 
