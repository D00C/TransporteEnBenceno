#!/bin/python
import TraF
import numpy as np
import csv
import matplotlib.pyplot as plt


def main():
    #Intervalo para la energía y número de iteraciones
    E = np.linspace(-1,1,1000)
    #Valores de la energía de sitio y los hoppings
    A = 0
    B = -0.5

    #Se abre el archivo para los datos
    outcsv = open("Transmission.csv","w")
    escritor = csv.writer(outcsv,delimiter=',',quotechar='',quoting=csv.QUOTE_NONE)
    while True:
        C = input('Dame la cadena de moleculas(p,m o q para salir): ')
        if C == 'q':
            break
        An,Bn = TraF.MOL(E,A,B,C)
        an,am,bnm = TraF.Nor(E,B,An,Bn,len(C))
        T = TraF.Tra(E,A,B,an,am,bnm)
        #for i in range(len(T)):
            #escritor.writerow([str(E[i]),str(T[i])])
    fig, ax = plt.subplots()
    ax.plot(E,T)
    plt.show()

if __name__ == '__main__':
    main()
