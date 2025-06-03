
from random import randint,random
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import lagrange
import struct 


def modifie_i_eme(message:str,i:int,new:str):
    """prend une chaine de caractère et renvoie la même chaine avec son ieme caractère changé en new

    Args:
        message (str): chaine
        i (int): indice
        new (str): nouveau caractère

    Returns:
        str: voir description
    """
    if i < 0 or i >= len(message):  # Vérifie que l'indice est valide
        return message
    return message[:i] + new + message[i+1:]

def repr_bin_str(c:str):
    """Convertit une chaîne en une représentation binaire


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

def uniforme(bit:str,proba:int):
    """procéde à l'opération binaire not sur le message avec une probabilité 1/p

    Args:
        bit (str): taille = 1
        proba (int): proba>0

    Returns:
        str: renvoie l'opposé du bits
    """
    p = randint(1, proba)
    if p ==1 : 
        if bit=='1':
            return '0'
        else:
            return '1'
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
    #offset_period=randint(0,period-1)
    assert loi in l 
    assert proba>1
    n = len(message)
    i=0
    while i < n : 
        if loi == 'uniforme' :
            message = modifie_i_eme(message,i,uniforme(message[i],proba))
        elif loi == 'periodique':
            assert period>0 
            if i%period == 0: 
                message = modifie_i_eme(message,i,uniforme(message[i],1))
        elif loi == 'burst_uni' :
            assert taille_burst>0 
            if 1== randint(1, proba) : 
                j = 0 
                while (i+j)<n and j<taille_burst :
                    message = modifie_i_eme(message,i,uniforme(message[i],1))
                    j = j +1 
                i=i+64
        i=i+1
    return message 

def est_puissance_de_2(n:int):
    """vérifie si n est une puissance de 2

    Args:
        n (int): un entier lambda

    Returns:
        bool: true ssi n est une puissance de 2
    """
    return n > 0 and (n & (n - 1)) == 0 # & et bit a bit

def decoupage_puissance_2(chaine:str,m:int):  
    """découpe un str en bloc de taille 2^m et rajoute des bits de parité pour correspondre au codage de hamming

    Args:
        chaine (str): str
        m (int): >=4

    Returns:
        list: liste de str de taille 2^m 
    """
    assert m>=3
    n = len(chaine)
    j = 0 
    p = 2 **m
    c = ''.join([' ' for i in range(p)])
    longueur=n//(p - m -1)
    if longueur*11<n : 
        longueur= longueur +1
    t=[c for i in range(longueur)]
    h = 0
    while j<n : 
        c = 'p'
        if j + p -1-m< n :
            for i in range(1,p):
                if est_puissance_de_2(i): # rajoute bits de redondance
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
            for l in range(i,p): # on rajoute des 0 à la fin pour avoir un message de taille une puissance de 2
                if est_puissance_de_2(l):
                    c= c + 'r'
                else:
                    c = c + '0'
                    j = j +1 
            if h<len(t) and len(t[h])<2**m : 
                t[h]=c+'0'
            else:
                t[h]=c 
            h=h+1
    return t 

     
def calculate_parity_positions(r:int):
    """renvoie les positions des bits de parité dans un code de hamming

    Args:
        r (int): entier

    Returns:
        list: liste contenant les positions
    """
    # les positions de la case de redondance i sont ceux telles que le i eme bits de la case représenté est un 1
    n = 2**r - 1  
    parity_positions = {i: [] for i in range(1, r+1)}
    for i in range(1, n + 1):
        for p in range(1, r + 1):
            if (i & (2**(p - 1))): 
                parity_positions[p].append(i)
    return parity_positions

def xor_bit(c1:str,c2:str):
    """xor

    Args:
        c1 (str): '0' ou '1'
        c2 (str): '0' ou '1'

    Returns:
        str: '0' ou '1'
    """
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
    """renvoie la somme (xor) des positions des bits '1' dans le message cela permet de trouver la valeur du bits de parité hamming

    Args:
        pos (list): liste de position
        message (str): chaine de '0' ou '1'

    Returns:
        str: '0' ou '1'
    """
    c = '0'
    for elt in pos:
        c = xor_bit(c,message[elt])
    return c

def redondance(message:str,m:int):
    """donne une valeurs aux bits de redondance

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
    """Bit de détection d'erreur Hamming (bit 2)

    Args:
        message (str): str
        m (int): pour faire une puissance de 2

    Returns:
        str: le message avec son premier caractère modifier pour hamming
    """
    n = 2**m
    compte=0
    for i in range(1,n):
        if message[i]==1:
            compte+=1
    if compte%2==0:
        return '0'+message[1:]
    else:
        return '1'+message[1:]

def position_erreur(message:str):
    """
    trouve la position de l'erreur s'il n'y en a qu'une
    sinon renvoie 0 
    Args:
        message (str): str
    Returns:
        int: indice de l'erreur
    """
    t=[ i for i,bit in enumerate(message) if bit=='1']
    pos = 0
    for elt in t : 
        pos =pos ^ elt 
    return pos 
    
def correction_erreur(message:str): 
    """corrige une erreur 

    Args:
        message (str): str

    Returns:
        str: message initial 
    """
    i = position_erreur(message)
    if i == 0 or i>=len(message):
        return message
    return modifie_i_eme(message,i,uniforme(message[i],1))

def hamming_encoding(message:str,m:int):
    """codage de hamming du message avec des decoupages de taille 2^m

    Args:
        message (str): message initial
        m (int): >=4

    Returns:
        list: liste de str de taille exactement 2^m représentant le message initial envoyé avec la redondance
    """
    c= codage_binaire(message) # transforme en binaire
    #c='10001111010'
    l = decoupage_puissance_2(c,m) # on découpe en plus petit message de taille 2**m
    n = len(l)
    for i in range(n):
        l[i] = redondance(l[i],m) # calcul des bits de redondance
        l[i] = detecteur_1(l[i],m) # calcul du bit de dectection
    return l
#print(position_erreur(hamming_encoding('',4)[0]))
def position_not_p2(m:int):
    """pour un entier m donné renvoie un tableau contenant tous les entiers de 1 à m qui ne sont pas des puissances de 2

    Args:
        m (int): >=1

    Returns:
        list: voir description
    """
    pos=[]
    for i in range(1,m):
        if not est_puissance_de_2(i):
            pos.append(i)
    return pos

def hamming_decoding(l:list):
    """decode un message qui a subit une transformation de hamming

    Args:
        l (list): list de str de même taille contenant uniquement des 0 ou 1 

    Returns:
        str: message decoder qui sera le message initial si il n'y a qu'une erreur 
    """
    message=''
    n= len(l)
    m=len(l[0])
    pos=position_not_p2(m)
    for i in range(n):
        l[i]=correction_erreur(l[i]) # corrige erreur de chaque morceau du message
        for indice in pos:
            message=message + l[i][indice] # on recupère le message sans les bits de redondance/parité
    p = len(message) 
    h=p//7 
    m=''
    for i in range(h): # on reconstruit le message ASCII
        c=''
        for j in range(7):
            if message[7*i+j]!=' ': # ?
                c=c+message[7*i+j]
            else:
                c=c+'0'
        if len(c)==7:
            asci_i = int(c,2)
            if asci_i!=0:
                m=m+ chr(asci_i)
    return m

#Reed-Solomon

def float64_to_bin(value:float):
    """Convertit un float64 en binaire IEEE 754 (64 bits)."""
    packed = struct.pack('!d', value)  # '!d' = Big-endian double (64 bits)
    return ''.join(f'{byte:08b}' for byte in packed)

def bin64_to_float(binary:str):
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
            valeur = float('inf')  
        else:
            valeur = float('nan')  
    else:
        # Nombre normalisé
        valeur = (1 + mantisse / (2**52)) * 2**(exposant - 1023)

    # Appliquer le signe
    return -valeur if signe else valeur

def valeur_str(chaine:str):
    """
    renvoie  pour chaine = m1m2m3m4 ...
    Somme des mi*2**i

    Args:
        chaine (str): chaine de 0 et 1

    Returns:
        int: 
    """
    v= 0 
    p = 1
    for elt in chaine:
        if elt=='1':
            v = v + p
        p=2*p 
    return v 

def inverse_valeur_str(v: int,p:int):
    """inverse de la fonction précédente

    Args:
        v (int): la valeur renvoyer par la dernière fonction
        p (int): définir la taille de la chaine renvoyer (rajoute des 0 non significatif)

    Returns:
        str: 
    """
    if v == 0:
        return "0".ljust(p, "0")  # Retourne une chaîne de 14 zéros si v == 0
    
    chaine = ""
    while v > 0:
        chaine += "1" if v % 2 == 1 else "0"
        v //= 2
    
    return chaine.ljust(p, "0")  # Complète avec des zéros à droite jusqu'à 14 caractères

def decoupage_taille_p(message:str,p:int):
    """découpe un str en bloc de taille p 

    Args:
        message (str): str
        p (int): >0 ici p = 64

    Returns:
        list: liste de str de taille p
    """
    n = len(message)
    q= n//p
    r = n%p
    l = ['' for i in range(p+1)]
    for i in range(p+1):
        c=''
        for j in range(q):
            if i*q + j < n : 
                c=c+ message[i*q + j]
        l[i]=c
    if r != 0 :
        l.append(message[(p)*q+q:])
    return l 

def encodage_reed_solomon(message:str,redondance:int,p:int):
    """
    Code de Reed-Solomon avec des paquets.

    Args:
        message (str): message quelconque.
        redondance (int): nombre de caractères de redondance.
        p (int): taille des paquets.

    Returns:
        list: liste des paquets.
    """
    message=codage_binaire(message) # codage binaire
    l = decoupage_taille_p(message,p) # découpage pour valeurs des polynomes
    n = len(l)
    m=n//2
    x_points=np.array([i for i in range(n)])
    for i in range(m,n):
        x_points[i]+=redondance # ajouts de la redondance
    y_points=np.array([np.float64(0) for i in range(n)])
    for i in range(n):
        y_points[i]=np.float64(valeur_str(l[i])) 
        l[i]=y_points[i]
    P=lagrange(x_points,y_points) # construction du polynome
    for r in range(redondance): # rajoute les valeurs de redondance
        l.append(P(m+r))#float64_to_bin(P(m+r)))
    return l 

def correction_perte(l,k):
    """retrouve le message originale avec Polynome d'interpolation 

    Args:
        l (list): de float64
        k (int): nombre de caractères de redondance

    Returns:
        list: messsage
    """
    n = len(l)
    m=(n-k)//2
    x_points=np.array([i for i in range(n-k)])
    for i in range(m,n-k):
        x_points[i]+=k
    y_values=['False' for i in range(n)]
    c=0
    for i in range(n):# on retrouve les valuers non perdues
        if l[i]!='':
            y_values[i]=(l[i]) 
        else:
            c=c+1
    if c>k:
        return []
    else:
        x=[]
        y=[]
        dec=0
        compteur=0
        while compteur<n-k: # construction des tableaux x et y pour retrouver le polynome
            if y_values[compteur]!='False':
                x.append(x_points[compteur])
                y.append(y_values[compteur])
            else:
                x.append(m+dec)#x_points[compteur]+m-compteur+dec
                y.append(y_values[n-k+dec])
                dec+=1
            compteur+=1
        if 'False' in y:
            return []
        P=lagrange(x,y)
        l_message=['' for i in range(n-k)]###on retrouve le message originale
        m=(n-k)//2
        i=0
        j=0
        while i<n:
            l_message[j]=(np.rint((P(x_points[j]))))
            if i+1 ==m:
                i=i+k+1
            else:
                i=i+1
            j=j+1
        return l_message 
    
def retrouvemessageascii(l,p):
    """reconstruction du message original

    Args:
        l (lst str): message bianire
        p (int): taille des blocs

    Returns:
        str: message
    """
    nn=len(l)
    n = p*nn
    m=n//7
    message_bin=''
    message=''
    ll=l.copy()
    for i in range(nn):
        ll[i]=inverse_valeur_str((ll[i]),p)
        message_bin+=ll[i]
    for i in range(m):
        mm=''
        if 7*i + 6 <n :
            for j in range(7):
                mm+=message_bin[7*i+j]
        message+=chr(int(mm,2))
    mm=''.join([elt for elt in message if elt!='\x00' and elt!='\n' and elt!='\x0c' and elt!='\x05' and elt!='\x0b'])
    return mm
            
def random_disapear(l,proba:int):
    """probabilité de perdre le contenu d'un elt de la liste

    Args:
        l (lst str): liste
        proba (int): plus grand que 1 

    Returns:
        lst str: potentiellement des messages devenues vide
    """
    n = len(l)
    for i in range(n) : 
        if random() < 1/proba : 
            l[i]=''
    return l

def test_RS(m:str,int:int,proba:int,redondance:int):
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
        n=len(codage_binaire(m))//15 + 1 
        l=encodage_reed_solomon(m,redondance,n)
        ll=random_disapear(l,proba) #probabilité d'enlever un bloc de 64 bits
        l2= correction_perte(ll,redondance)
        mm=retrouvemessageascii(l2,n)
        if mm==m:
            c=c+1
        if i%(int/10)==0 and False:
            print(i)
    return c 

def test_sans_cor(m:str,int:int,proba:int):
    """

    Args:
        m (str): message
        int (int): nombre de fois que l'on test le procedé
        taille (int): sur combien de bits utilisé hamming
        proba (int): 1/proba = probabilité qu'un bits soit alteré (voir dans la fonction si uniforme ou pas)
    
    Returns:
        int: nombre de transmission correcte 
    """
    c=0
    for i in range(int+1):
        mm= random_change(m,'uniforme',proba)
        if mm==m:
            c=c+1
    return c 

def test_hamming(m:str,int:int,taille:int,proba:int,perte:int):
    """Fonction de test du codage de Hamming dans un canal bruité.
    On teste `int` fois la transmission et on compte le nombre de fois où le message est correctement retranscrit.


    Args:
        m (str): message
        int (int): nombre de fois que l'on test le procedé
        taille (int): sur combien de bits utilisé hamming
        proba (int): 1/proba = probabilité qu'un bits soit alteré (voir dans la fonction si uniforme ou pas)
        perte(int): = 0 corruption de donnée = 1 perte de donnée

    Returns:
        int: nombre de transmission correcte 
    """
    c=0
    corrompu= ''.join(['0' for i in range(2**taille)])
    for i in range(int):
        l=hamming_encoding(m,taille)
        n=len(l)
        for i in range(n):
            if perte == 1 and random()<1/proba: # dans le cas d'une perte de donné 
                l[i]= corrompu
            elif perte==0:
                l[i]= random_change(l[i],'uniforme',proba)
        mm=hamming_decoding(l)
        if mm==m:
            c=c+1
    return c

def entrelacement(l: list, p: int = 16):
    """
    Entrelace les bits de la liste `l` par blocs de `p` bits.
    Chaque bloc est vu comme une ligne, et l'entrelacement se fait colonne par colonne.
    """
    n = len(l)
    m = n // p  # nombre de blocs
    # Construire la matrice m x p (m lignes de p bits)
    # Entrelacement : lire colonne par colonne
    new_l = []
    for k in range(m):
        for i in range(p): # pour chaque colonne
            mess=''
            for j in range(p):  # pour chaque ligne
                mess=mess+l[p*k+j][i]
            new_l.append(mess)
            
    for i in range(p*m,p*m+n%p):
        new_l.append(l[i])
    return new_l

message='This is a simple ASCII text '
n=len(codage_binaire(message))//15 + 1 
#print(len(hamming_encoding(message,6))* len(hamming_encoding(message,6)[0]))
#print(len(encodage_reed_solomon(message,4,n))* 64)
l1=[''.join([chr(j) for i in range(16)]) for j in range(97,97+16)]
l2=[''.join(format(j,'016b') ) for j in range(16)]
l1= l1 + l1

def perte_random(l):
    """
    change un elt de la liste en 16 '0'
    Args:
        l (list): str

    Returns:
        list: str 
    """
    n = len(l)
    j=randint(1,n-1)
    l[j]=''.join(['0' for i in range(len(l[j]))])
    return l

m='This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text aaaaa'

assert entrelacement(entrelacement(l1))==l1
assert m == hamming_decoding(entrelacement(entrelacement(hamming_encoding(m,4))))
assert m == hamming_decoding(entrelacement(perte_random(entrelacement(hamming_encoding(m,4))))) # si 16 | len(hamming_encoding(m,4))


def poly_add(p1:list, p2:list, modulo:int):
    """ addition de 2 polynomes modulo"""
    length = max(len(p1), len(p2))
    result = [0] * length
    for i in range(length):
        if i < len(p1):
            a=p1[i]
        else :
            a=0
        if i < len(p2):
            b=p2[i]
        else :
            b=0
        result[i] = (a + b) % modulo
    return result

def poly_mul(p1:list, p2:list, modulo:int):
    # multiplication de deux polynomes
    result = [0] * (len(p1) + len(p2) - 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            result[i + j] = (result[i + j] + p1[i] * p2[j]) % modulo
    return result

def lagrange_poly(x:list, y:list, modulo=131):
    # renvoie sous forme d'une liste de coeffs le polynome d'interpolation de lagrange modulo modulo
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


def format_polynomial(coeffs:list,modulo:int):
    #pour afficher le polynôme
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


def eval_poly(coeffs:list, x:int, modulo:int):
    # renvoie P(x) où P est décrit par sa liste de ses coeeficients
    result = 0
    for coef in reversed(coeffs):
        result = (result * x + coef) % modulo
    return result



def encodage_reed_solomon_fini(message:str,redondance:int,modulo=131):
    """
    Code de Reed-Solomon avec des paquets avec corps fini

    Args:
        message (str): message quelconque.
        redondance (int): nombre de caractères de redondance.
        modulo (int): modulo

    Returns:
        list: liste des paquets.
    """
    n = len(message)
    l = decoupage_taille_p(message,n)
    ll=[ 0 for i in range(n+redondance)]
    x_points=[i for i in range(n)] # i+1 pour le test de la source
    y_points=[ord(message[i]) for i in range(n)]
    for i in range(n):
        ll [i]=y_points[i]

    P=lagrange_poly(x_points,y_points,modulo)
    for r in range(redondance):
        ll[n+r]=eval_poly(P,n+r,modulo)
    return ll 



def correction_perte_fini(l,k,modulo=131):
    """retrouve le message originale avec Polynome d'interpolation modulo modulo

    Args:
        l (list): message potentiellement corrompu
        k (int): nombre de caractères de redondance

    Returns:
        list: messsage
    """
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
    #renvoie le codage ASCII d'une liste d'entiers
    return ''.join([chr(elt) for elt in l])


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
    plt.ylabel("Nombre de messages correctement retrouvés")
    plt.grid()
    plt.show()



# test 2**m avec m qui change pour hamming
def hamming_pour_different_n():
    m='Hello! This is a simple ASCII text with numbers 1234567890 and symbols: @#$%^&*()-_=+[]{};:,.<>?'
    print(len(codage_binaire(m)))
    l=hamming_encoding(m,3)
    print(len((l)))
    x_values=[4,5,6,7,8]
    y_values=[]
    m='Hello! This is a simple ASCII text with numbers'
    m=m+'1234567890 and symbols: @#$%^&*()-_=+[]{};:,.<>?'
    mp = test_sans_cor(m,1000,500)
    yy_values=[mp for i in range(len(x_values))]
    for i in range(len(x_values)):
        y_values.append(test_hamming(m,1000,x_values[i],500,0))
        print(i)
    plt.scatter(x_values,y_values) # points bleu
    plt.xlabel(" valeur de n")
    plt.scatter(x_values,yy_values) #points orange
    plt.ylabel("Nombre de messages correctement retrouvés")
    plt.grid()
    plt.show()

#hamming_pour_different_n()

#  test niveau taille avec n=4 

# test entrelacement Hamming vs RS 

def test_hamming_entrelace(m:str,int:int,taille:int,proba:int):
    """fonction de test du codage de hamming dans un canal bruité
    on test int fois la fonction et on compte le nombre de fois où le message est correctement retranscrit 

    Args:
        m (str): message
        int (int): nombre de fois que l'on test le procedé
        taille (int): sur combien de bits utilisé hamming
        proba (int): 1/proba = probabilité qu'un bits soit alteré (voir dans la fonction si uniforme ou pas)
        perte(int): = 0 corruption de donnée = 1 perte de donnée

    Returns:
        int: nombre de transmission correcte 
    """
    c=0
    corrompu= ''.join(['0' for i in range(2**taille)]) # si il y avait un 0 pas d'erreur et sinon il la dectetera
    for i in range(int):
        l=hamming_encoding(m,taille)
        n=len(l)
        l = entrelacement(l,2**taille)
        for i in range(n):
            if random()<1/proba: # dans le cas d'une perte de donné 
                l[i]= corrompu
        l = entrelacement(l,2**taille)
        mm=hamming_decoding(l)
        if mm==m:
            c=c+1
    return c


#Hamming_vs_Reed_Solomon()
# test valeur de RS


def Reed_Solomon_different_k():
    message='This is a simple ASCII text '
    m='This is a simple ASCII text '
    n=len(codage_binaire(m))//15 + 1 
    ll=encodage_reed_solomon(message,3,n)
    print(len(ll))
    print('test')
    x_values=[]
    y_values=[]
    nb=1000
    m='This is a simple ASCII text '
    for i in range(0,4):
        print(i)
        x_values.append(i)
        y_values.append(test_RS(message,nb,50,i))
    print(sum(y_values)/len(y_values))
    plt.scatter(x_values,y_values)
    plt.xlabel("valeur de k")
    plt.ylabel("Nombre de messages correctement retrouvés")
    plt.grid()
    plt.show()
    
#Reed_Solomon_different_k()


def Hamming_vs_Reed_Solomon():
    message_h='This is a simple ASCII text '
    message_RS='This is a simple ASCII text '
    n=len(codage_binaire(message_RS))//15 + 1 
    l=hamming_encoding(message_h,6)
    ll=encodage_reed_solomon(message_RS,3,n)
    print(len(l))
    print(len(ll))
    nb=10000
    redondance_RS=6
    taille_hamming=6
    
    nh=test_hamming_entrelace(message_h,nb,taille_hamming,20)
    nrs=test_RS(message_RS,nb,20,redondance_RS)
    categories = ['Hamming' , 'Reed-Solomon']
    values= [nh,nrs]
    colors=['blue','red']
    print(values)
    plt.bar(categories,values,color=colors)
    plt.ylim(0,nb)
    plt.xlabel("")
    plt.ylabel("Nombre de messages correctement retrouvés")
    plt.grid()
    plt.show()

        
Hamming_vs_Reed_Solomon()

def Hamming_vs_Reed_Solomon_fini():
    message_h='This is a simple ASCII text '
    message_RS='This is a simple ASCII text '
    l=hamming_encoding(message_h,3)
    ll=encodage_reed_solomon_fini(message_RS,4)
    print(l)
    print(len(l))
    print(ll)
    print(len(ll))
    nb=1000
    redondance_RS=5
    taille_hamming=3
    
    nh=test_hamming_entrelace(message_h,nb,taille_hamming,35.555)
    nrs=test_RS_fini(message_RS,nb,35.555,redondance_RS)
    categories = ['Hamming' , 'Reed-Solomon']
    values= [nh,nrs]
    colors=['blue','red']
    print(values)
    plt.bar(categories,values,color=colors)
    plt.ylim(0,nb)
    plt.xlabel("")
    plt.ylabel("Nombre de messages correctement retrouvés")
    plt.grid()
    plt.show()
    
#Hamming_vs_Reed_Solomon_fini()