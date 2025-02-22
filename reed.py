"""#Reed-Solomon

import numpy as np
from scipy.interpolate import lagrange

# Points d'interpolation (x, y)
x_points = np.array([1, 2, 3])
y_points = np.array([2, 3, 5])

# Calcul du polynôme interpolateur
#P = lagrange(x_points, y_points)

#print(P)  # Affichage du polynôme

def valeur_str(chaine:str):
    v= 0 
    p = 1
    for elt in chaine:
        if elt=='1':
            v = v + p
        p=2*p 
    return v 

def decoupage_taille_p(message:str,p:int):
    l = ['' for i in range(p)]
    n = len(message)
    m= n//p
    for i in range(p):
        c=''
        for j in range(m):
            c=c+ message[i*p + j]
        l[i]=c
    if n %p != 0 :
        l.append(message[(p-1)*p+m:])
    return l 

def encodage_reed_solomon(message:str,redondance:int,p:int):
    message=codage_binaire(message)
    l = decoupage_taille_p(message,p)
    n = len(l)
    x_points=np.array([i for i in range(n)])
    y_points=np.array([i for i in range(n)])
    for i in range(n):
        y_points[i]=valeur_str(l[i])

    print(y_points)
    P=lagrange(x_points,y_points)
    for r in range(redondance):
        print(P(n+r))
        #message+= bin(P(n+r))[2:] # forcer le 7 bits
        l.append(bin(int(P(n+r)))[3:])
    return l 

print(encodage_reed_solomon('elisa',3,5))
"""