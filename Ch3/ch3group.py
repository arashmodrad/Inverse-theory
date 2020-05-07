import numpy as np
from scipy.io import loadmat
from matplotlib import pyplot as plt

colscan = loadmat("colscan.mat")['colscan']
rowscan = loadmat("rowscan.mat")['rowscan']
diag1 = loadmat("diag1scan.mat")['diag1scan']
diag2 = loadmat("diag2scan.mat")['diag2scan']

# form the g3 matrix 
n = 16
m = 31 # diag m dimension = 2n-1


# right diagonal matrix; diag1 
def rdg(m,n):
    G3 = np.zeros((m,n*n))
    G3[0,0] = 1             # its super annoying to do this in the loop below 
    for i in range(n):
        G3[i, list(np.arange(0,n*i,n-1) + i)] = 1 
        if i<n-1:
            G3[i+n, list(np.arange(n*n-(n-i-1),n*(i+1),-n+1))] = 1 
        """
        Rows are shaped like ... where x denotes number of 1s in the column 
             x 
             x x 
             x x x
             x x x x
             x x x
             x x
             x 
        """
    return G3

# left diagonal matrix; diag2 
def ldg(m,n):
    # left diagonal matrix 
    # diag 2
    n = 16
    m = 31
    G4 = np.zeros((m,n*n))
    n=16
    for i in range(n):
        G4[i, list(np.arange(n-1+n*i, -1, -(n+1)))] = 1
        if i < n-1:
            G4[i+n, list(np.arange(n**2-i-2, n*i, -(n+1)))] = 1 
    # 
    # 
    return G4


#------------------------------# 
# create matrices 
G1 = np.zeros((n, n**2)) # rowscan
for i in range(16):
    G1[i, i*n:(i+1)*n] = 1 


G2 = np.zeros((n,n**2))    # colscan
for i in range(16):
    G2[i, list(np.arange(i, n**2, n))] = 1

G3 = rdg(m,n)*np.sqrt(2) # diag1

G4 = ldg(m,n)*np.sqrt(2) # diag2 

# Construct full matrix 
def stack_matrices(mlist):
    # mlist is a list of matrices
    from functools import reduce
    rows = reduce(lambda x,y: x+y, [i.shape[0] for i in mlist])
    cols = [i.shape[1] for i in mlist]
    if cols.count(cols[0]) != len(cols):
        print "cols do not match size"
        return 
    # matrix to output 
    stack = np.zeros((rows, cols[0]))
    lr=0
    for mat in mlist:
        r=mat.shape[0]+lr
        stack[lr:r,:] = mat
        lr = r
    return stack 

def stack_vectors(vlist):
    rows = reduce(lambda x,y: x+y, [i.shape[0] for i in vlist])
    stack = np.zeros((rows, 1))
    lr = 0
    for vec in vlist:
        r=vec.shape[0]+lr
        stack[lr:r,0] = vec
        lr = r 
    return stack 
# rowscan
# colscan
# diag1
# diag2


# First part --- just row/col vectors
Ga = stack_matrices([G1, G2])
t1 = stack_vectors([rowscan[:,-1], colscan[:,-1]])

#u,s,v = np.linalg.svd(G1)
#p = s.size - s[s<1e-10].size
#print 'n,m', G.shape, 'p=', p 

# all of the data  
Gb = stack_matrices([G1, G2, G3, G4])
t2 = stack_vectors([rowscan[:,-1], colscan[:,-1], diag1[:,-1], diag2[:,-1]])
#u,s,v = np.linalg.svd(G)
#p = s.size - s[s<1e-10].size
#print 'n,m', G.shape, 'p=', p 



# calculate the model resoultion matrices
def mRes(G):
    u,s,v = np.linalg.svd(G)
    v = v.T  # fix V... ugh why do they do this 
    # how many ps 
    p = s.size - s[s<1e-10].size 
    
    # select p columns
    vp = v[:,0:p]
    
    # model res matrix 
    R = vp.dot(vp.T)
    # daig only
    diagR = np.diagonal(R)
    # reshape
    resz = int(diagR.size**.5)

    rshp = np.reshape(diagR, (resz, resz))
    
    return rshp

# make a plot 
def mResplot():
    fig, ax= plt.subplots(1,2)
    i1 = ax[0].imshow(mRes(Ga), vmin=.2, vmax=1.5)
    ax[0].set_title("Matrix A")

    i2 = ax[1].imshow(mRes(Gb), vmin=.2, vmax=1.5)
    ax[1].set_title("Matrix B")

# find the pseudoinvere solution
from ind3 import pinv
mdag1 = pinv(Ga, 31).dot(t1).reshape((16,16))
mdag2 = pinv(Gb, 87).dot(t2).reshape((16,16))

def nullplot(G1,G2,plot):
    u1,s1,v1 = np.linalg.svd(G1)
    u2,s2,v2 = np.linalg.svd(G2)
    v1 = v1.T
    v2 = v2.T
    p1 = s1.size - s1[s1<1e-10].size
    p2 = s2.size - s2[s2<1e-10].size
    # create nullspace plots 
    if plot=="yes":
        plt.figure(1)
        plt.subplot(121)
        plt.imshow(v1[:,p2+1].reshape((16,16)))
        plt.title("Model Null Space for A")
        
        plt.subplot(122)
        plt.imshow(mdag2)
        plt.imshow(v2[:,p2+1].reshape((16,16)))
        plt.title("Model Null Space for B")
        plt.colorbar()
    return v1[:,p1+1].reshape((16,16)), v2[:,p2+1].reshape((16,16))



def slowplot():
    plt.figure(1)
    plt.subplot(121)
    plt.imshow(mdag1)
    plt.title("min: {} \n max: {}".format(np.round(mdag1.min(),8), np.round(mdag1.max(),8)))

    
    plt.subplot(122)
    plt.imshow(mdag2)
    plt.title("min: {} \n max: {}".format(np.round(mdag2.min(),8), np.round(mdag2.max(),8)))
    plt.show()

#v1,v2 = nullplot(Ga,Gb,'no')
def wildplot():
    plt.figure(1)
    plt.subplot(121)
    plt.imshow(mdag1+v1*20, vmin=-5, vmax=20)
#    plt.title("min: {} \n max: {}".format(np.round(mdag1.min(),8), np.round(mdag1.max(),8)))
    
    plt.subplot(122)
    plt.imshow(mdag2+v2*20.,vmin=-5, vmax=20)
 #   plt.title("min: {} \n max: {}".format(np.round(mdag2.min(),8), np.round(mdag2.max(),8)))
    plt.show()



