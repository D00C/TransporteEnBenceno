def dim2(E,A,B):
    #Constantes para simplificar la ecuación
    X = [E-A[0],E-A[1],E-A[2],E-A[3]]
    z = X[1]*X[2] - B[1]**2
    #
    A1 = A[0] + X[2]*B[0]**2/z
    A4 = A[3] + X[1]*B[2]**2/z
    B14 = B[0]*B[1]*B[2]/z
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
        an,am,bnm = dim2(E,An[-1][0:4],Bn[-1][0:3])
        An.append([an,am])
        Bn.append([bnm])
        i = 1
        while i <= (n//2 - 1):
            #Se agrega el hopping de la unión
            Bn[-1].append(Bn[-2][4*i-1])
            #Se toman los dimeros n y n+1 y se normalizan a uno solo en el ultimo indice del arreglo
            an,am,bnm = dim2(E,An[-2][4*i:4*(i+1)+1],Bn[-2][4*i:4*(i+1)])
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

def MOL(E,a,b,Cad,D=None):
    X = (E-a)/b
    #Para-benceno
    Bp = 2*b/(X**2-1)
    Ap = a + Bp*X
    #Meta-benceno
    Bm = b*(X**2-1)/(X*(X**2-2))
    Am = a + b/X + Bm
    #Orto-benceno
    Bo = b + b/(X**2*(X**2 - 2))
    Ao1 = a + b*(X**2 - 1)/(X*(X**2 - 2))
    Ao2 = a + b*(X**6 - 2*X**4 + X**2 - 1)/(X**3*(X**2 - 1)*(X**2 - 2))
    A = [[]]
    B = [[]]
    if D == None:
        D = []
        for i in range(len(Cad) - 1):
            D.append(1)
    for i in range(len(Cad)):
        if Cad[i] == 'p':
            A[0].append(Ap)
            A[0].append(Ap)
            B[0].append(Bp)
        elif Cad[i] == 'm':
            A[0].append(Am)
            A[0].append(Am)
            B[0].append(Bm)
        elif Cad[i] == 'o':
            A[0].append(Ao1)
            A[0].append(Ao2)
            B[0].append(Bo)
        else:
            print('Error, simbolo incorrecto: '+i)
            break
        if i < len(Cad)-1:
            B[0].append(D[i]*b)
    return A,B
