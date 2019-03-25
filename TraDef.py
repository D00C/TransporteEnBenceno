#!/bin/python
import TraF
import numpy as np
import pandas as pd

def main():
    #Intervalo de deformación y número de iteraciones
    D = np.linspace(0,0.1,1000)
    #D = np.linspace(1,1,1000)
    #Intervalo para la energía y número de iteraciones
    E = np.linspace(-1,1,1000)
    #Valores de la energía de sitio y los hoppings
    A = 0
    B = -0.5

    #String que define las configuraciones de la cadena
    #p es configuracion para
    #m es configuracion meta
    #o es configuracion ortho
    C = "mmmmmm"
    d = []
    for i in range(len(C) - 1):
        d.append(1)
    T = []
    for i in range(len(D)):
        for j in range(len(d)):
            if j < len(C)/2:
                d[j] = d[j] + D[i]*((j+1)/len(C))
            else:
                d[j] = d[j] + D[i]*(1 - (j+1)/len(C))
        #Se crea los arreglos de energias de sitio y de enlace
        An,Bn = TraF.MOL(E,A,B,C,d)
        #Se renormaliza la molecula
        an,am,bnm = TraF.Nor(E,B,An,Bn,len(C))
        #Se calcula la transmicion
        T.append(TraF.Tra(E,A,B,an,am,bnm))
    df = pd.DataFrame(data=T,index=D,columns=E)
    df.to_csv(path_or_buf='TD.csv')

if __name__ == '__main__':
    main()
    import plot
