import matplotlib.pyplot as plt
from math import *
import numpy as np
import csv

def DIM2(E,B,a,b):
    #Constantes para simplificar la ecuación
    x = [(E-a[0])/B,(E-a[1])/B,(E-a[2])/B,(E-a[3])/B]
    N = [b[0]/B,b[1]/B,b[2]/B]
    z = x[1] - N[1]**2/x[2]

    #Formulas para la normalización
    #A1 = a1 + N1²B/z
    A1 = a[0] + B*N[0]**2/z
    #A4 = a4 + N3²B/x3 + (N2N3/X3)²B/z
    A4 = a[3] + B*N[2]**2/x[2] + (N[1]*N[2]/x[2])**2*B/z
    #B14 = N1N2N3B/(X3z)
    B14 = B*N[0]*N[1]*N[2]/(x[2]*z)
    return A1,A4,B14

def Tra(E,A,B,An,Am,Bnm):
    #Constantes para simplificar la ecuación
    zn = (An-A)/(2*B)
    zm = (Am-A)/(2*B)
    G = (Bnm-B)/(2*B)
    P = zn + zm
    Q = zn*zm - G - G**2
    X = (E-A)/B

    #Formula para la transmisión
    #T = (1+2G)²(4-X²)/[(1-2Q)²(4-X²) + 4(P-QX)²]
    T = (1+2*G)**2*(4-X**2)
    T /= (1-2*Q)**2*(4-X**2)+4*(P-Q*X)**2
    return T

def Nor(E,B,An,Bn,n):
    while n > 1:
        #Se toman los primeros dos dimeros y se normalizan a uno solo en un nuevo indice del arreglo
        an,am,bnm = DIM2(E,B,An[-1][0:4],Bn[-1][0:3])
        An.append([an,am])
        Bn.append([bnm])
        i = 1
        while i <= (n//2 - 1):
            #Se agrega el hopping de la unión
            Bn[-1].append(Bn[-2][4*i-1])
            #Se toman los dimeros n y n+1 y se normalizan a uno solo en el ultimo indice del arreglo
            an,am,bnm = DIM2(E,B,An[-2][4*i:4*(i+1)+1],Bn[-2][4*i:4*(i+1)])
            An[-1].append(an)
            An[-1].append(am)
            Bn[-1].append(bnm)
            i += 1
        #En caso de que quede un ultimo dimero se coloca al final del ultimo indice del arreglo
        if n%2 != 0:
            An[-1].append(An[-2][-2])
            An[-1].append(An[-2][-1])
            Bn[-1].append(Bn[-2][-2])
            Bn[-1].append(Bn[-2][-1])
        n -= n//2
    return An[-1][-2],An[-1][-1],Bn[-1][-1]

def Mol(Amol,Bmol,B,n):
    if isinstance(Amol,list):
        pass
    else:
        #Se crea el arreglo con las 2n energías de sitio
        An = [[]]
        for i in range(2*n):
            An[0].append(Amol)

        #Determina los coeficientas de deformacion para el cable
        #if n > 1:
            #Ans = input('¿Hay deformaciones?(s/n): ')
        #else:
            #Ans = None
        #k = []
        #if Ans == 's':
            #for i in range(n-1):
                #ki = input("Dame el coeficiente k%i : "%(i+1))
                #try:
                    #k.append(float(ki))
                #except (ValueError):
                    #k.append(eval(ki))
        #else:
            #for i in range(n-1):
                #k.append(1)

        #Se crea arreglo con las 2n-1 energías de enlace
        Bn = [[]]
        for i in range(n-1):
            Bn[0].append(Bmol)
            #Bn[0].append(B*k[i])
            Bn[0].append(B)
        Bn[0].append(Bmol)
    return An,Bn


def main():
    #Intervalo para la energía y número de iteraciones
    E = np.linspace(-1,1,8000)
    #Valores de la energía de sitio y los hoppings
    A = 0
    B = -0.5

    X = (E-A)/B
    #Para-benceno
    Bp = 2*B/(X**2-1)
    Ap = A + Bp*X
    #Meta-benceno
    Bm = B*(X**2-1)/(X*(X**2-2))
    Am = A + B/X + Bm

    while True:
        fig, ax = plt.subplots()
        n = int(input('Número de moleculas: '))
        An,Bn = Mol(Am,Bm,B,n)
        an,am,bnm = Nor(E,B,An,Bn,n)
        T = Tra(E,A,B,an,am,bnm)
        nom = input('Nombre del perfil: ')

        ###Gráfica
        #Reducir el tamaño del plot un 20%
        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #Colocar la leyenda a la derecha del plot
        #ax.legend(bbox_to_anchor=(1,0.5))
        ax.plot(E,T)
        ax.grid()
        plt.title(nom+' moleculas')
        plt.xlabel('E')
        plt.ylabel('T')
        plt.show()
        Ans = input("¿Parar el programa?(s/n): ")
        if Ans == 's':
            break

if __name__ == '__main__':
    main()
