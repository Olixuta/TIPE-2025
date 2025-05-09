import numpy as np
from scipy.interpolate import lagrange
x=[i for i in range(10)]
print(x)
y=[1429,13259,9602.5,359,1,8,3,0,-5,9.0]
P=(lagrange(x,y))
x=[i for i in range(1,11)]
print(x)
y=[13259,9602.5,359,1,8,3,0,-5,9.0,P(10)]
P2=lagrange(x,y)

print(P)
print(P2)

X1=[1,2,3,4,5]
Y1=[7,8,-3,7,17]
P1=lagrange(X1,Y1)
X2=[2,3,4,5,6]
Y2=[8,-3,7,17,P1(6)]
P2=lagrange(X2,Y2)
print(P1)



# decodage reed solomon
def random_disapear(l,proba:int):
    n = len(l)
    for i in range(n) : 
        if random() < 1/proba : 
            l[i]=''
    return l
print(P2)