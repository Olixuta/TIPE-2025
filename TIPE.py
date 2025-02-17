 

# essayer vite fait de faire une modélisation de hamming pour étudier les différents type d'erreurs qui pourrait advenir
# périodique / en groupe / uniforme 


# on envoie les bits d'un même mot dans un ordre différents avec les autres mots pour éviter des erreurs par paquets. 




# par exemple codons huffman pour avoir des codes de base sur une chaine


#faire une fonction tri

"""
def huffman (chaine:str):
    n = len(chaine) 
    d= {}
    for elt in chaine :
        if not elt in d : 
            d[elt] = 1 
        else:
            d[elt] = d[elt] + 1 
    print(dict.items(d)) 
    for key in d: 
        print(key)
        
        

        
huffman('bbvc')"""

import numpy as np 
from functools import reduce
import operator as op
from random import randint
import matplotlib as plt


def modifie_i_eme(message:str,i:int,new:str):
        return message[:i] + new + message[i+1:]



def repr_bin_str(c:str):
    """converti un charactère en sa représentation binaire

    Args:
        c (str): un str de taille 1 

    Returns:
        str: représentation binaire de c 
    """
    return bin(ord(c))[2:]


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

def uniforme(bit:str,proba:int):
    """procéde à l'opération binaire not sur le message avec une probabilité 1/p

    Args:
        bit (str): taille = 1
        proba (int): proba>0

    Returns:
        str: description
    """
    p = randint(1, proba)
    if p ==1 : 
        print('m')
        if bit=='0':
            return '1' 
        else:
            return '0'
    return bit



def random_change(message:str,loi:str,proba:int,period:int=0,taille_burst:int=0):
    """modélise un canal bruité 

    Args:
        message (str): message initial
        loi (str): type d'erreurs expected: 'uniforme','periodique','burst_uni'
        proba (int): probabilité qu'un bit ou qu'un ensemble de bits soit altéré 
        period (int, optional): si les erreurs arrivent de manière périodique. Defaults to 0.
        taille_burst (int, optional): taille d'un bloc. Defaults to 0.

    Returns:
        str: le message avec probablement des erreurs 
    """
    l=['uniforme','periodique','burst_uni']
    assert loi in l 
    assert proba>1
    n = len(message)
    i=0
    while i < n : 
        if loi == 'uniforme' : 
            message = modifie_i_eme(message,i,uniforme(message[i],proba))
        elif loi == 'periodique':
            assert period>0 
            if i%period == 0 : 
                message = modifie_i_eme(message,i,uniforme(message[i],1))
        elif loi == 'burst_uni' :
            assert taille_burst>0 
            if 1== randint(1, proba) : 
                j = 0 
                while (i+j)<n and j<taille_burst :
                    message = modifie_i_eme(message,i,uniforme(message[i],1))
                    j = j +1 
        i=i+1
    return message 



def est_puissance_de_2(n):
    return n > 0 and (n & (n - 1)) == 0

def decoupage_puissance_2(chaine:str,m:int): 
    n = len(chaine)
    j = 0 
    p = 2 **m
    c = ''.join([' ' for i in range(p)])
    t=[c for i in range(n//(p - m -1) + 1)]
    h = 0
    while j<n : 
        c = 'p'
        if j + p -1-m< n :
            for i in range(1,p):
                if est_puissance_de_2(i):
                    c=c+ 'r'
                else:
                    c =c+ chaine[j]
                    j = j +1 
            t[h]=c
            h=h+1
        else: 
            i = 1 
            while j <n:
                if est_puissance_de_2(i):
                    c= c + 'r'
                else:
                    c = c + chaine[j]
                    j = j +1 
                i = i +1
            for l in range(i,p):
                if est_puissance_de_2(l):
                    c= c + 'r'
                else:
                    c = c + '0'
                    j = j +1 
            t[h]=c+'0'
            h=h+1
    return t 

     

def calculate_parity_positions(r):
    # les positions de la case de redondance i sont ceux telles que le i eme bits de la case représenté est un 1
    n = 2**r - 1  # Total number of bits in Hamming code (m + r)
    # A dictionary to store which bits each parity covers
    parity_positions = {i: [] for i in range(1, r+1)}
    for i in range(1, n + 1):
        # For each bit position i, check which parity bits cover it
        for p in range(1, r + 1):
            if (i & (2**(p - 1))):  # Check if bit p covers position i
                parity_positions[p].append(i)
    return parity_positions


def xor_bit(c1:str,c2:str):# je sais pas si c'est un xor
    if c1=='1':
        if c2=='1':
            return '0'
        else:
            return '1'
    else:
        if c2=='1':
            return '1'
        else:
            return '0'



def parity_pos(pos:list,message:str):
    c = '0'
    for elt in pos:
        c = xor_bit(c,message[elt])
    return c


def redondance(message:str,m:int):
    """choisit les bits de redondance

    Args:
        message (str): chaine de taille m
        m (int): entier 
    """
    d = calculate_parity_positions(m)
    for key in d: 
        indice = 2**(key-1 )
        message=modifie_i_eme(message,indice,parity_pos(d[key],message))
    return message
        
def detecteur_1(message:str,m:int):
    c = '0'
    n = 2**m
    for i in range(1,n):
        c = xor_bit(c,message[i])
    return c+message[1:]



def position_erreur(message:str):
    t=[ i for i,bit in enumerate(message) if bit=='1']
    return reduce(op.xor , t)
    
def correction_erreur(message:str): 
    i = position_erreur(message)
    if i == 0:
        return message
    return modifie_i_eme(message,i,uniforme(message[i],1))



#je t\'aime
#tu as reussi ton ds de maths ?
#j\'ai envie de t\'embrasser
#elisa
#love
#Hello! This is a simple ASCII text with numbers 1234567890 and symbols: @#$%^&*()-_=+[]{};:,.<>?

def hamming_encoding(message:str,m:int):
    c= codage_binaire(message)
    l = decoupage_puissance_2(c,m)
    n = len(l)
    for i in range(n):
        l[i] = redondance(l[i],m)
        l[i] = detecteur_1(l[i],m)
    return l 


l=hamming_encoding('elisa',4)
print(l)

l[0]= random_change(l[0],'uniforme',50)
print(l)
print(position_erreur(l[0]))