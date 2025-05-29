#Reed-Solomon

import numpy as np
from scipy.interpolate import lagrange
import struct
import numpy as np
from random import random

    

def float64_to_bin(value):
    """Convertit un float64 en binaire IEEE 754 (64 bits)."""
    packed = struct.pack('!d', value)  # '!d' = Big-endian double (64 bits)
    return ''.join(f'{byte:08b}' for byte in packed)

def bin64_to_float(binary):
    """Convertit une chaîne binaire IEEE 754 (64 bits) en float64."""
    assert len(binary) == 64, "La chaîne doit contenir exactement 64 bits"

    # Extraction des parties
    signe = int(binary[0], 2)
    exposant = int(binary[1:12], 2)
    mantisse = int(binary[12:], 2)

    # Calcul du float64
    if exposant == 0:
        # Nombre subnormal (exposant = 0 → exposant réel = -1022)
        valeur = (mantisse / (2**52)) * 2**-1022
    elif exposant == 2047:
        # Cas spéciaux : Infini ou NaN
        if mantisse == 0:
            valeur = float('inf')  # Infini
        else:
            valeur = float('nan')  # NaN
    else:
        # Nombre normalisé
        valeur = (1 + mantisse / (2**52)) * 2**(exposant - 1023)

    # Appliquer le signe
    return -valeur if signe else valeur


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


def valeur_str(chaine:str):
    v= 0 
    p = 1
    for elt in chaine:
        if elt=='1':
            v = v + p
        p=2*p 
    return v 

def inverse_valeur_str(v: int,p:int) -> str:
    if v == 0:
        return "0".ljust(p, "0")  # Retourne une chaîne de 14 zéros si v == 0
    
    chaine = ""
    while v > 0:
        chaine += "1" if v % 2 == 1 else "0"
        v //= 2
    
    return chaine.ljust(p, "0")  # Complète avec des zéros à droite jusqu'à 14 caractères




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


def encodage_reed_solomon(message:str,redondance:int,p:int):
    message=codage_binaire(message)
    l = decoupage_taille_p(message,p)
    n = len(l)
    m=n//2
    x_points=np.array([i for i in range(n)])
    for i in range(m,n):
        x_points[i]+=redondance
    y_points=np.array([np.float64(0) for i in range(n)])
    for i in range(n):
        #print(valeur_str(l[i]))
        y_points[i]=np.float64(valeur_str(l[i]))
        l[i]=y_points[i]
    #print(x_points)
    P=lagrange(x_points,y_points)
    '''
    for i in range(n):
        l[i]=y_points[i]#float64_to_bin(valeur_str(l[i]))
        #print(type(y_points[i]))
        #print(type(l[i]))'''
    for r in range(redondance):
        #print(P(m+r))
        #assert(abs(P(m+r))<=2**15)
        l.append(P(m+r))#float64_to_bin(P(m+r)))
    #print(x_points)
    #print(len(l))
    #print(P(7))
    #print(P)
    return l 





def correction_perte(l,k):#plus que cette fonction
    n = len(l)
    m=(n-k)//2
    x_points=np.array([i for i in range(n-k)])
    for i in range(m,n-k):
        x_points[i]+=k
    print(x_points)
    x_values=[i for i in range(n)]
    #print(x_values)
    #print(x_points)
    y_values=['False' for i in range(n)]
    c=0
    for i in range(n):# on retrouve les valuers non perdues
        if l[i]!='':
            #print(type(l[i]))
            y_values[i]=(l[i]) #np.float64(valeur_str(l[i]))
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
                print('   ')
                print(dec)
                print(compteur)
                print(n-k+dec)
                x.append(x_points[compteur]+m-compteur+dec) # -1 car l[1]='' etc voir pouquoi sauf que ça marche plus avec plusieurs erreurs
                y.append(y_values[n-k+dec])
                dec+=1
            compteur+=1
        #print(x)
        #print(y)
        print(x)
        print(y)
        if 'False' in y:
            return [] # pas possible
        P=lagrange(x,y)###
        #print(P)
        l_message=['' for i in range(n-k)]###on retrouve le message originale
        m=(n-k)//2
        i=0
        j=0
        print('t')
        print(x_points)
        print(x)
        while i<n:
            #print(x_points)
            #print(j)
            #print(x_points[j])
            l_message[j]=(np.rint((P(x_points[j]))))#normalement np.rint
            #"print(type(l_message[i]))# chr(valeur_str(float64_to_bin
            if i+1 ==m:
                i=i+k+1
            else:
                i=i+1
            j=j+1
            #print(j)
        ###
        return l_message 
    
def retrouvemessageascii(l,p):
    print(l)
    nn=len(l)
    n = p*nn
    m=n//7
    message_bin=''
    message=''
    ll=l.copy()
    for i in range(nn):
        #print(type(l[i]))
        ll[i]=inverse_valeur_str((ll[i]),p)
        message_bin+=ll[i]
    for i in range(m):
        mm=''
        if 7*i + 6 <n :
            for j in range(7):
                mm+=message_bin[7*i+j]
        #print(mm)
        #print(int(mm,2))
        message+=chr(int(mm,2))
    mm=''.join([elt for elt in message if elt!='\x00' and elt!='\n' and elt!='\x0c' and elt!='\x05' and elt!='\x0b'])
    return mm
        
            
            
            
def random_disapear(l,proba:int):
    n = len(l)
    for i in range(n) : 
        if random() < 1/proba : 
            l[i]=''
    return l


def test_RS(m:str,int:int,proba:int,redondance:int):
    """_summary_

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
        n=len(codage_binaire(m))//15 + 1 
        l=encodage_reed_solomon(m,redondance,n)
        ll=random_disapear(l,proba) #probabilité d'enlever un bloc de 64 bits
        l2= correction_perte(ll,redondance)
        mm=retrouvemessageascii(l2,n)
        if mm==m:
            c=c+1
        #print(i)
        if i%100==0:
            print(i)
    return c 

redondance=3

m='This is a simple ASCII text '
n=len(codage_binaire(m))//15 +1 
l1=(encodage_reed_solomon(m,redondance,n)) # rajoute 3 paquets 
l1[3]='' # perte de la partie 1 et 2 du message
l1[2]=''
l1[1]='' # 
#l1[2]=''
l2=(correction_perte(l1,redondance))
message =retrouvemessageascii(l2,n)
assert message == m
n=len(codage_binaire(m))//15 +1 # en modifiant la taille de m ou cette ligne ça marche pas

#print(test_RS(m,100,10,3))
#print(test_RS(m,100,1000,3))
#print(l1)
#print(l2)