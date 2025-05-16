
import operator as op
from random import randint,random
from matplotlib import pyplot as plt
from time import time 
import numpy as np
from scipy.interpolate import lagrange
import struct
from math import ceil,floor  


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
        #print('m')
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

def decoupage_puissance_2(chaine:str,m:int):  # à vérifier pour ' '
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
    #print(len(t))
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
m= codage_binaire("This is a simple ASCII code")
#print(m)
#print(len(m))
#print(decoupage_puissance_2(m,4))
     

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
    """bit 2 detecteur d'erreur hamming

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
    """trouve la position de l'erreur si il n'y en a qu'une
    sinon renvoie 0 (ce bit est utilisé pour 2 détecteur)

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
    #print(l[0])
    n = len(l)
    for i in range(n):
        l[i] = redondance(l[i],m) # calcul des bits de redondance
        l[i] = detecteur_1(l[i],m) # calcul du bit de dectection
    return l

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
    message=codage_binaire(message) # codage binaire
    l = decoupage_taille_p(message,p) # découpage pour valeurs des polynomes
    n = len(l)
    m=n//2
    x_points=np.array([i for i in range(n)])
    for i in range(m,n):
        x_points[i]+=redondance # ajouts de la redondance
    y_points=np.array([np.float64(0) for i in range(n)])
    for i in range(n):
        #print(valeur_str(l[i]))
        y_points[i]=np.float64(valeur_str(l[i])) 
        l[i]=y_points[i]
    #print(x_points)
    P=lagrange(x_points,y_points) # construction du polynome
    for r in range(redondance): # rajoute les valeurs de redondance
        #assert(abs(P(m+r))<=2**64)
        l.append(P(m+r))#float64_to_bin(P(m+r)))
    return l 



def correction_perte(l,k):
    n = len(l)
    m=(n-k)//2
    x_points=np.array([i for i in range(n-k)])
    for i in range(m,n-k):
        x_points[i]+=k
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
        return 'rip'
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
                x.append(x_points[compteur]+m)
                y.append(y_values[n-k+dec])
                dec+=1
            compteur+=1
        #print(x)
        #print(y)
        P=lagrange(x,y)###
        #print(P)
        l_message=['' for i in range(n-k)]###on retrouve le message originale
        m=(n-k)//2
        i=0
        j=0
        while i<n:
            l_message[j]=(np.rint((P(x_points[j]))))
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
    """probabilité de perdre le contenu d'un elt de la liste

    Args:
        l (lst str): _description_
        proba (int): _description_

    Returns:
        lst str: _description_
    """
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
        if i%(int/10)==0:
            print(i)
    return c 

m='This is a simple ASCII text '
'''l1=(encodage_reed_solomon(m,3,n))
print(l1)
print(retrouvemessageascii(l1,n))
l2=(correction_perte(l1,3))
print(l2)
print(retrouvemessageascii(l2,n))'''
n=len(codage_binaire(m))//15 +1 # en modifiant la taille de m ou cette ligne ça marche pas
#print(test_RS(m,1000,1000,3))



def test_sans_cor(m:str,int:int,proba:int):
    """_summary_

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
    corrompu= ''.join(['0' for i in range(2**taille)])
    for i in range(int+1):
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

def entrelacement(l:list,p=16):
    #la liste sort de hamming encoding
    # entrelace p blocs de p bits (rajouter en paramètre pour modifié)
    n = len(l)
    m = n //p
    print('im here')
    print(m)
    l_matrice=[[['' for i in range(p)]for k in range(p)] for j in range(1)] # m+n%p
    dec = 0 
    for k in range(p):
        for j in range(m):
            for i in range(p):
                l_matrice[j][k][i]+= l[p*j+i][k]
    print(l_matrice)
    new_l=['' for i in range(n)]
    for j in range(m):
        for i in range(p):
            for k in range(p):
                new_l[j+k]+= l_matrice[j][k][i]
    return new_l 

message='This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text '
l1=[''.join([chr(j) for i in range(16)]) for j in range(97,97+16)]
l2=[''.join(format(j,'016b') ) for j in range(16)]
l1.append('a')
print(l1)
print(entrelacement(l1))
#print(entrelacement(entrelacement(l1)))
#print((hamming_encoding('10001111010',4)[0]))
#01010000111110100


# une chance sur 1000 d'erreur uniforme avec bloc de taille 16 -> 9985
#(testsurtaphrasemignongne) 16 -> 9976 32 -> 9944 64->9891 128->9770 256->9391 512->9071
# proba 100 -> (512 -> 360(36secondes) 256(32secondes),16->7637(21secondes) )
#t1=time()
#print(len(m))
#print(test_hamming(m,10000,4,100))
#print(time()-t1)
# test niveau temps

'''
# test 2**m avec m qui change pour hamming
m='Hello! This is a simple ASCII text with numbers 1234567890 and symbols: @#$%^&*()-_=+[]{};:,.<>?'
print(len(codage_binaire(m)))
l=hamming_encoding(m,3)
print(len((l)))
x_values=[3,4,5,6,7,8,9,10]
y_values=[]
mp = test_sans_cor(m,1000,500)
print(mp)
yy_values=[mp for i in range(len(x_values))]
time_values=[]
t1=time()
for i in range(len(x_values)):
    y_values.append(test_hamming(m,1000,x_values[i],500,0))
    t2=time()
    time_values.append(100*(t2-t1))
    t1=t2
    print(i)
plt.scatter(x_values,y_values) # bleu
plt.scatter(x_values,yy_values) # bleu
plt.show()'''


"""
# test m=4 et on augmente la taille du message
m='Hello! This is a simple ASCII text with numbers 1234567890 and symbols: @#$%^&*()-_=+[]{};:,.<>?'
print(test_hamming(m,100,4,100,0))
x_values=[]
y_values=[]
yy_values=[]
yyy_values=[]
rh=[]
rsc=[]
y_sanc_cor=0
y_h=0
t1=time()
for i in range(0,1000,20):
    y_h+=test_hamming(m,20,3,100,0)
    y_sanc_cor+=test_sans_cor(m,20,100)
    x_values.append(i)
    y_values.append(y_h)
    yy_values.append(y_sanc_cor)
    yyy_values.append(i)
    t2=time()
    t1=t2
    #print(i)
plt.scatter(x_values,y_values) # bleu
plt.scatter(x_values,yy_values)
plt.scatter(x_values,yyy_values)
print(len(codage_binaire(m)))
print(y_values[-1]/1000)
print(yy_values[-1]/1000)
plt.show()"""

#  test niveau taille avec n=4 

"""# test corruption de donnée (marche pas)
print(len(codage_binaire(m)))
message='This is a simple ASCII text '
x_values=[i for i in range(4)]
y_values=[]
yy_values=[]
yyy_values=[]
time_values=[]
t1=time()
t0=time()
for i in range(len(x_values)):
    y_values.append(test_hamming(message,100*i,8,400)) # -> une chance sur 100 pour un bloc de taille 16
    yy_values.append(test_sans_cor(message,100*i,100))
    yyy_values.append(test_RS(message,100*i,156,3))
    t2=time()
    time_values.append(100*(t2-t1))
    t1=t2
    print(i)
    
    
print(time()-t0)
plt.scatter(x_values,yy_values) # bleu
plt.scatter(x_values,y_values) #orange
plt.scatter(x_values,yyy_values) #vert
plt.scatter(x_values,[100*i for i in range(len(x_values))]) # rouge
#plt.scatter(x_values,time_values)
plt.show()"""

""" Test perte de donné 
message='This is a simple ASCII text '
print(len(codage_binaire(message)))
n=len(codage_binaire(m))//15 + 1 
l=hamming_encoding(message,6)
ll=encodage_reed_solomon(message,3,n)
print(len(l))
print(len(ll))
print('test')
print(1000-test_hamming(message,1000,6,100,1))#151
print(1000-test_RS(message,1000,100*9/2,3))#19 block de taille plus faible donc plus de proba d'alterner moduler pour cette cause
""" 



""" les tailles des burst modifiées peuvvent variée de qlq bits à qlq centaines"""
#fibres qq dizaines maxes




# proba que sur des blocs de 16 avec une proba d'erreurs de 1/100 on est 2 erreurs : 0.0109
# proba que sur des blocs de 32 avec une proba d'erreurs de 1/200 on est 2 erreurs : 0.0112
# proba que sur des blocs de 64 avec une proba d'erreurs de 1/400 on est 2 erreurs : 0.0114



#reed solomon même proba 

# test entrelacement Hamming vs RS 

def test_hamming_entrelace(m:str,int:int,taille:int,proba:int): # marche pas vérifie
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
    print(len(corrompu))
    for i in range(int+1):
        l=hamming_encoding(m,taille)
        n=len(l)
        print(l)
        l = entrelacement(l,2**taille)
        print(l)
        for i in range(n):
            if random()<1/proba: # dans le cas d'une perte de donné 
                l[i]= corrompu
        print(l)
        l = entrelacement(l,2**taille)
        mm=hamming_decoding(l)
        if mm==m:
            c=c+1
    return c


message='This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text This is a simple ASCII text '
print(len(decoupage_puissance_2(message,4)))
print(test_hamming_entrelace(message,100,4,10)) 
# plt.grid() -> avant le plt.show()