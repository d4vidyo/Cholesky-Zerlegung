import numpy as np

def toLaTeX(A):
    #Konvertiere zu LaTeX
    lines = str(A)[1:-1].replace('[','').replace('\n','').split(']')
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)

class Band:
    def __init__(self, A, m, symetrie = "symetrisch"):
        n = len(A)
        self.m = m
        self.symetrie = symetrie

        self.Bmatrix = np.zeros((m+1, n), dtype = A.dtype)
        for i in range(0, m+1):
            self.Bmatrix[i][i:n] = np.diagonal(A, i)
            

    def getElem(self, zeile, spalte):
        if self.symetrie == "symetrisch":
            if zeile > spalte: return self.getElem(spalte, zeile)
            index = spalte - zeile
            if index > self.m: return 0
            return self.Bmatrix[index][spalte]
        elif self.symetrie == "unteresDreieck":
            if zeile < spalte: return 0
            index = zeile - spalte
            if index > self.m: return 0
            return self.Bmatrix[index][zeile]

    def setElem(self, zeile, spalte, value):
        if zeile > spalte: return self.setElem(spalte, zeile, value)
        index = spalte - zeile
        if index > self.m: return
        self.Bmatrix[index][spalte] = value
        
    def ausschreiben(self):
        n = len(self.Bmatrix[0])
        A = np.zeros((n,n), dtype = self.Bmatrix.dtype)
        for i in range(0, n):
            for j in range(0, n):
                A[i][j] = self.getElem(i,j)
        return A

    def erzeugeZufaellige(n = None):
        n = n if n else np.random.randint(low = 8, high = 16)
        m = np.random.randint(low = 1, high = n-1)

        Rand = np.zeros((n,n))
        randDiag = np.random.random(n) * 10
        Rand += np.diag(randDiag, k=0)
        for i in range(1,m+1):
            randDiag = np.random.random(n-i) - 0.5 * np.ones(n-i)
            Rand += np.diag(randDiag,k=i)
            
        A = Rand.T.dot(Rand)

        if not np.all(np.linalg.eigvalsh(A)>0): return Band.erzeugeZufaellige()
        
        return Band(A, m)


    def __str__(self):
        return str(self.Bmatrix)

    def __getattr__(self, attr):
        return getattr(self.Bmatrix, attr)
    

def choleskyComp(A, m, **options):
    typ = options.get('dtype') or np.float64
    n = len(A.Bmatrix[0])
    L = Band(np.zeros((n,n), dtype = typ), m , "unteresDreieck")
    for i in range(0,n):
        Sum = typ(0)
        for k in range(0,i):       
            Sum += np.square(L.getElem(i, k))      
        L.setElem(i, i, np.sqrt(A.getElem(i, i) - Sum))        

        for j in range(i+1, n):
            Sum = typ(0)
            for k in range(0, i):
                Sum += L.getElem(j, k) * L.getElem(i, k)
            L.setElem(j, i, (A.getElem(j, i) - Sum) / L.getElem(i, i))

    return L

def main():
    typ = np.float64
    epsilon = 0.0001

    zufallstest = True
    for i in range(0,999):
        A = Band.erzeugeZufaellige()
        L = choleskyComp(A, A.m)

        L_lang = L.ausschreiben()        
        Rekonstruiert = L_lang.dot(L_lang.T)
        A_lang = A.ausschreiben()
        if not np.all(A_lang - Rekonstruiert < epsilon):
            zufallstest = False
            break

    print("Zufallstest: ", zufallstest)


    B = np.array([
        [4,-1,0,0],
        [-1,4,-1,0],
        [0,-1,4,-1],
        [0,0,-1,4]
    ], dtype = typ)
    I = np.eye(4, dtype = typ)
    Z = np.zeros((4,4), dtype = typ)

    A1 = np.hstack((B, -I, Z, Z))
    A2 = np.hstack((-I, B, -I, Z))
    A3 = np.hstack((Z, -I, B, -I))
    A4 = np.hstack((Z, Z, -I, B))
    A = np.vstack((A1, A2, A3, A4))

    m = 4
    n = len(A)
    
    A_kompakt = Band(A, m)

    L = choleskyComp(A_kompakt, m, dtype = typ)
    
    L_lang = L.ausschreiben()
    Rekonstruiert = L_lang.dot(L_lang.T)
    print("Vorlesungsbeispiel Test: ", np.all(A - Rekonstruiert < epsilon))#
    print("L=" + toLaTeX(np.round(L,1)) + "\\\\")
    print("L_{lang}="+toLaTeX(np.round(L_lang,1)))


if __name__ == "__main__":
    main()    
