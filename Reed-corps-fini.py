#Reed-Solomon

import operator as op
from random import randint,random,choice
from matplotlib import pyplot as plt
from time import time 
import numpy as np
from scipy.interpolate import lagrange
import struct 


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

m= 'This is a simple ASCII code '
n=len(m)
print(decoupage_taille_p(m,int(n/2)))
    
    
    

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

def lagrange_poly(x, y, modulo=131):
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


def eval_poly(coeffs, x, modulo):
    result = 0
    for coef in reversed(coeffs):
        result = (result * x + coef) % modulo
    return result



def encodage_reed_solomon_fini(message:str,redondance:int,modulo=131):# ok
    n = len(message)
    l = decoupage_taille_p(message,n)
    ll=[ 0 for i in range(n+redondance)]
    x_points=[i for i in range(n)] # i+1 pour le test de la source
    y_points=[ord(message[i]) for i in range(n)]
    #print('le message ascii')
    #print(y_points)
    for i in range(n):
        ll [i]=y_points[i]

    P=lagrange_poly(x_points,y_points,modulo)
    #print(P)
    for r in range(redondance):
        #print(P(m+r))
        #assert(abs(P(m+r))<=2**15)
        ll[n+r]=eval_poly(P,n+r,modulo)

    return ll 



def correction_perte_fini(l,k,modulo=131):
    n = len(l)
    x_points=([i for i in range(n)])
    y_values=['False' for i in range(n)]
    c=0
    for i in range(n):# on retrouve les valuers non perdues
        if l[i]!='':
            y_values[i]=(l[i]) 
        else:
            c=c+1
    if c>k:
        return [] # pas possible
    else:
        x=[]
        y=[]
        dec=0
        compteur=0
        while compteur<n-k:
            if y_values[compteur]!='False':
                x.append(x_points[compteur])
                y.append(y_values[compteur])
            else:
                x.append(x_points[n-k+dec])
                y.append(y_values[n-k+dec])
                dec+=1
            compteur+=1
        if 'False' in y:
            return [] # pas possible
        P=lagrange_poly(x,y,modulo)
        l_message=['' for i in range(n-k)]
        j=0
        while j<n-k:
            l_message[j]=eval_poly(P,x_points[j],modulo)
            j=j+1
        return l_message 
    
def decodage_ascii_fini(l):
    m = ''
    for elt in l:
        m+=chr(elt)
    return m


ll=encodage_reed_solomon_fini('This is a simple ASCII code ',4,131)
print(len(ll))
print(ll)
print(ll)
ll[0]=''
ll[1]=''
ll[2]=''
l3= correction_perte_fini(ll,4,131)
print(l3)
print(decodage_ascii_fini(l3))

def random_disapear(l,proba:int):
    """probabilité de perdre le contenu d'un elt de la liste
    """
    n = len(l)
    for i in range(n) : 
        if random() < 1/proba : 
            l[i]=''
    return l

def test_RS_fini(m:str,int:int,proba:int,redondance:int):
    """
    Args:
        m (str): message
        int (int): nombre de fois que l'on test le procedé
        taille (int): sur combien de bits utilisé hamming
        proba (int): 1/proba = probabilité qu'un bits soit alteré (voir dans la fonction si uniforme ou pas)
        redondance(int): nombre de redondance pour le polynome
    
    Returns:
        int: nombre de transmission correcte 
    """
    c=0
    for i in range(int):
        l=encodage_reed_solomon_fini(m,redondance)
        ll=random_disapear(l,proba) #probabilité d'enlever un bloc 
        l2= correction_perte_fini(ll,redondance)
        mm=decodage_ascii_fini(l2)
        if mm==m:
            c=c+1
        if i%(int/10)==0 and False:
            print(i)
    return c 

def Reed_Solomon_fini_different_k():
    message='This is a simple ASCII text '
    #print(len(codage_binaire(message)))
    ll=encodage_reed_solomon_fini(message,3)
    print(len(ll))
    print('test')
    x_values=[]
    y_values=[]
    nb=1000
    for i in range(0,7):
        print(i)
        x_values.append(i)
        y_values.append(test_RS_fini(message,nb,10,i))
    plt.scatter(x_values,y_values)
    plt.xlabel("valeur de k")
    plt.ylabel(" nombre de message correctement retrouvé")
    plt.grid()
    plt.show()


#Reed_Solomon_different_k()

