

# essayer vite fait de faire une modélisation de hamming pour étudier les différents type d'erreurs qui pourrait advenir
# périodique / en groupe / uniforme 


# on envoie les bits d'un même mot dans un ordre différents avec les autres mots pour éviter des erreurs par paquets. 




# par exemple codons huffman pour avoir des codes de base sur une chaine


#faire une fonction tri


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
        
        
        
huffman('bbvc')