
from functools import reduce
import operator as op
from random import randint
from matplotlib import pyplot as plt
from time import time 

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
        i=i+1
    return message 



def est_puissance_de_2(n:int):
    """vérifie si n est une puissance de 2

    Args:
        n (int): un entier lambda

    Returns:
        bool: true ssi n est une puissance de 2
    """
    return n > 0 and (n & (n - 1)) == 0

def decoupage_puissance_2(chaine:str,m:int): 
    """découpe un str en bloc de taille 2^m et rajoute des bits de parité pour correspondre au codage de hamming

    Args:
        chaine (str): str
        m (int): >=4

    Returns:
        list: liste de str de taille 2^m 
    """
    assert m>=4
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

    Args:
        message (str): str

    Returns:
        int: indice de l'erreur
    """
    t=[ i for i,bit in enumerate(message) if bit=='1']
    return reduce(op.xor , t) if t else 0 
    
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

def position_not_p2(m:int):
    pos=[]
    for i in range(1,m):
        if not est_puissance_de_2(i):
            pos.append(i)
    return pos


def hamming_decoding(l:list):
    message=''
    n= len(l)
    m=len(l[0])
    pos=position_not_p2(m)
    for i in range(n):
        l[i]=correction_erreur(l[i])
        for indice in pos:
            message=message + l[i][indice]
    p = len(message) 
    h=p//7 # on enleve les bits de 0 qu'on a rajouter à la fin 
    m=''
    for i in range(h):
        c=''
        for j in range(7):
            if message[7*i+j]!=' ': # je comprends pas pourquoi ça arrive mais ça n'impacte pas les performances donc why not 
                c=c+message[7*i+j]
            else:
                c=c+'0'
        if len(c)==7:
            asci_i = int(c,2)
            if asci_i!=0:
                m=m+ chr(asci_i)
    return m


    

#m='Elisa je t aime de tout mon coeur tu me manques ma reine' 
m='Elisa je t aime de tout mon coeur tu me manques ma reine'

def test_hamming(m:str,int:int,taille:int,proba:int):
    c=0
    for i in range(int):
        l=hamming_encoding(m,taille)
        n=len(l)
        for i in range(n):
            l[i]= random_change(l[i],'uniforme',proba)
        mm=hamming_decoding(l)
        if mm==m:
            c=c+1
    return c

# une chance sur 1000 d'erreur uniforme avec bloc de taille 16 -> 9985
#elisajetaimedetoutmoncoeurtumemanquesmareinee 16 -> 9976 32 -> 9944 64->9891 128->9770 256->9391 512->9071
# proba 100 -> (512 -> 360(36secondes) 256(32secondes),16->7637(21secondes) )
#t1=time()
#print(len(m))
#print(test_hamming(m,10000,4,100))
#print(time()-t1)
mm='Hello! This is a simple ASCII text with numbers 1234567890 and symbols: @#$%^&*()-_=+[]{};:,.<>? do you like me ? my love i miss u how can i forget those night ? '
x_values=[4,5,6,7,8,9,10,11]
#y_values=[test_hamming(m,1000,x_values[i],100)for i in range(len(x_values))]
yy_values=[]
time_values=[]
t1=time()
for i in range(len(x_values)):
    yy_values.append(test_hamming(mm,1000,x_values[i],1000))
    t2=time()
    time_values.append(100*(t2-t1))
    t1=t2


plt.scatter(x_values,yy_values)
plt.scatter(x_values,time_values)
plt.show()













