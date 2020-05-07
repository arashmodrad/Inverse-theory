import numpy as np
from scipy.io import loadmat
from matplotlib import pyplot as plt 
from ind3 import pinv 

d = loadmat("ifk.mat")['d']
m = 20
n = m 
interval = (1,0)
G = np.zeros((m,n))

dx = 1./n

b = np.ones(20)
np.fill_diagonal(G, b)
np.fill_diagonal(G[:,1:], b[:-1])


def gxi(x,y):
    return x*np.exp(-x*y)

xi = np.arange(dx/2, 1, dx)

mesh = np.meshgrid(xi,xi)

G = gxi(mesh[0], mesh[1])

# create picard ratio plot
u,s,v = np.linalg.svd(G)
v = v.T
piclist = []
for i in range(n):
    pic = u[:,i].T.dot(d)/s[i]
    piclist.append(pic)

#tsvd solution
Gtsvd = pinv(G,8)
mdag = Gtsvd.dot(d)

# least squares solution
mglstq = np.linalg.lstsq(G,d)


res_mdag  = d - G.dot(mdag)
res_lstsq = d - G.dot(mglstq[0])
plt.plot(res_mdag)# 'mdag residual')
plt.plot(res_lstsq)# 'lstsq residual')

