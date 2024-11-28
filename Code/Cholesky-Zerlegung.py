import numpy as np

def bmatrix(a):
    #Konvertiere zu LaTeX
    lines = str(a)[1:-1].replace('[','').replace('\n','').split(']')
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)

class Band:
    def __init__(self, Matrix, m, symetrie = "symetrisch"):
        n = len(Matrix)
        self.m = m
        self.symetrie = symetrie

        self.Bmatrix = np.zeros((m+1, n), dtype = Matrix.dtype)
        for i in range(0, m+1):
            self.Bmatrix[i][i:n] = np.diagonal(Matrix, i)

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

    def getZeile(self, index):
        if index < 0: return np.zeros(n, dtype = self.Bmatrix.dtype)
        n = len(self.Bmatrix[0])
        p = np.zeros(n, dtype = self.Bmatrix.dtype)
        if self.symetrie == "symetrisch":
            p[index:index+self.m+1] = np.diagonal(self.Bmatrix, index)
            for i in range(index-1,min(0,index-1-self.m),-1):
                p[i] = self.getElem(index, i)
            return p
        elif self.symetrie == "unteresDreieck":
            for i in range(index-1,min(0,index-1-self.m),-1):
                p[i] = self.getElem(index, i)
            return p


    def setZeile(self, index, value):
        if index < 0: return
        n = len(self.Bmatrix[0])
        for i in range(0, n):
            self.setElem(index, i, value[i])
    
    def getSpalte(self, index):
        return self.getZeile(index)

    def setSpalte(self, index, value):
        self.setZeile(index, value)

    def __str__(self):
        return str(self.Bmatrix)

    def __getattr__(self, attr):
        return getattr(self.Bmatrix, attr)



def getElem(M, zeile, spalte):
    if zeile > spalte: return getElem(M, spalte, zeile)
    m = len(M)-1
    index = spalte - zeile
    if index > m: return 0
    return M[index][spalte]

def setElem(M, zeile, spalte, value):
    if zeile > spalte: return setElem(M, spalte, zeile, value)
    m = len(M)-1
    index = spalte - zeile
    if index > m: return 0
    M[index][spalte] = value

"""
def getZeile(M, index):
    if index < 0: return np.zeros(n, dtype = M.dtype)
    n = len(M[0])
    m = len(M)-1
    p = np.zeros(n, dtype = M.dtype)
    p[index:index+m+1] = np.diagonal(M, index)
    for i in range(index-1,min(0,index-1-m),-1):
        p[i] = getElem(M, index, i)
    return p


def setZeile(M, index, value):
    if index < 0: return
    n = len(M[0])
    m = len(M)-1
    for i in range(0, n):
        setElem(M, index, i, value[i])
    
def getSpalte(M, index):
    return getZeile(M, index)

def setSpalte(M, index, value):
    return setZeile(M, index, value)
"""

def cholesky(A, **options):
    typ = options.get('dtype') or np.float64
    n = len(A)
    L = np.zeros((n,n), dtype = typ)
    for i in range(0,n):
        Sum = typ(0)
        for k in range(0,i):
            Sum += np.square(L[i][k])            
        L[i][i] = np.sqrt(A[i][i] - Sum)        

        for j in range(i+1, n):
            Sum = typ(0)
            for k in range(0, i):
                Sum += L[j][k] * L[i][k]
            L[j][i] = (A[j][i] - Sum) / L[i][i]
            
    return L

def choleskyComp(A, m, **options):
    typ = options.get('dtype') or np.float64
    n = len(A.Bmatrix[0])
    X = Band(np.zeros((n,n), dtype = typ), m , "unteresDreieck")
    L = np.zeros((n,n), dtype = typ)
    for i in range(0,n):
        Sum = typ(0)
        Sum2 = typ(0)
        for k in range(0,i):
            Sum += np.square(L[i][k])            
            Sum2 += np.square(X.getElem(i, k))            
        L[i][i] = np.sqrt(A.getElem(i, i) - Sum)
        X.setElem(i, i, np.sqrt(A.getElem(i, i) - Sum2))        

        for j in range(i+1, n):
            Sum = typ(0)
            Sum2 = typ(0)
            for k in range(0, i):
                Sum += L[j][k] * L[i][k]
                Sum2 += X.getElem(j, k) * X.getElem(i, k)
            L[j][i] = (A.getElem(j, i) - Sum) / L[i][i]
            X.setElem(j, i, (A.getElem(j, i) - Sum2) / X.getElem(i, i))
            
    L2 = np.zeros((n,n), dtype = typ)
    for i in range(0,n):
        for j in range(0,n):
            L2[i][j] = X.getElem(i,j)
            

    print("Zerlegung compakt:")
    print(np.round(X))
    print("")
    print("Zerlegung entpackt:")
    print(np.round(L2))
    print("")
    print("Zerlegung normal:")
    print(np.round(L))
    print(np.all(L==L2))
    return L

def main():
    typ = np.float64    
    
    #A = np.array([
    #    [4,2,4,4],
    #    [2,10,5,2],
    #    [4,5,9,6],
    #    [4,2,6,9]
    #    ], dtype=typ)

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
    print("Matrix normal:")
    print(A)
    print("")


    X = np.zeros((m+1,n),dtype=typ)
    for i in range(0, m+1):
        X[i][i:n]=np.diagonal(A, i)
    X = Band(A, m)
    print("Matrix compakt:")
    print(X)
    print("")
    
    L = cholesky(A, dtype = typ)
    Lcomp = choleskyComp(X, m, dtype = typ)
    print(np.all(L==Lcomp))
    #print(bmatrix(np.round(L,1)) + "\\\\")
    #print(bmatrix(np.round(L2,1)))

    #print(bmatrix(np.round(L, 1)) + "\\\\")
    #print(bmatrix(np.round(Lcomp, 1)))
    #print(np.all(L == Lcomp))
    #print(np.round(L, 0))
    #print(np.round(L.dot(L.T),2))
    #print("")
    #print(np.all(abs(L.dot(L.T) - A) < 0.00001))
    #print(bmatrix(A))

if __name__ == "__main__":
    typ = np.float64    
    off_diag = [1, 2, 3]
    B = np.array([
        [4,-1,0,0],
        [-1,4,4,0],
        [0,-1,4,-1],
        [0,0,-1,4]
    ], dtype = typ)
    target_diag = B.diagonal(1)
    res_diag = off_diag - target_diag 
    B += np.diag(res_diag, k=1)
    main()    
