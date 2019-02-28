import TraF
import numpy as np
import csv


def main():
    #Intervalo de deformación y número de iteraciones
    D = np.linspace(-1,1,4)
    #Intervalo para la energía y número de iteraciones
    E = np.linspace(-1,1,1000)
    #Valores de la energía de sitio y los hoppings
    A = 0
    B = -0.5

    #Se abre el archivo para los datos
    outcsv = open("Transmission.csv","w")
    escritor = csv.writer(outcsv,delimiter=',',quotechar='',quoting=csv.QUOTE_NONE)
    #String que define las configuraciones de la cadena
    #p es configuracion para
    #m es configuracion meta
    #o es configuracion ortho
    C = "pppp"
    for i in range(len(D)):
        #Se crea los arreglos de energias de sitio y de enlace
        An,Bn = TraF.MOL(E,A,B,C,D[i])
        #Se renormaliza la molecula
        an,am,bnm = TraF.Nor(E,B,An,Bn,len(C))
        #Se calcula la transmicion
        T = TraF.Tra(E,A,B,an,am,bnm)
        for j in range(len(T)):
            escritor.writerow([str(E[j]),str(D[i]),str(T[j])])

if __name__ == '__main__':
    main()
