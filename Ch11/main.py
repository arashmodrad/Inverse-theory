import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# read data
data = pd.read_csv("soils_data.csv")
h = np.array(data["preshead"]).reshape(data.preshead.shape[0],1)
#theta = np.array(data["theta"]).reshape(data.theta.shape[0],1)

def VanGenuchten(h,x):   # functools partial takes the 1st fx argument
    # x0 = theta_r
    # x1 = theta_s
    # x2 = alpha
    # x3= n
   # x4 = m
    bottom = (1. + (x[2]*h)**x[3])**(x[4])
    phi = x[0] + (x[1] - x[0])/bottom
    return phi

m_real = np.array([.92, .4, .5, 1.2, 2.1])
theta = np.array([VanGenuchten(hval, m_real) for hval in list(h)])

def G(x,**kwargs):
#    h = kwargs.get("h", np.linspace(0,20.,80)) # default h data
    phi = np.array([VanGenuchten(hbar,x) for hbar in list(h)])
    return phi

class uniform():
    #
    def __init__(self, val, vmin,vmax):
        self.val  = val
        self.vmin = vmin
        self.vmax = vmax
        return

def logprior(mlist):
    lpl = []
    for m in mlist:
        if (m.val >= m.vmin) and (m.val <= m.vmax):
            lpl.append(0)
        else:
            lpl.append(-1*np.inf)
    # return 1 value; if any params are out of the range, return -inf
    if lpl.count(0) != len(lpl):
        return -1*np.inf
    else:
        return 0 

def logliklihood(G, m, data, sigma):
    #
    fvec = (data - G(m, h=h))/sigma
    l = (-1/2)*np.sum(fvec**2)
    return l 

def logproposal(x,y,step):
    return (-1/2)*np.sum((x-y)**2/step**2)

def Generate(m0,step):
    return m0 +  step*np.random.randn(m0.shape[0], m0.shape[1])

    
def MC(**kwargs):
    vmin = 0  
    vmax = 1.   # max val for the parameters
    m0 = np.random.rand(5,1)
#def MCMC(m0):
    niter = kwargs.get("niter", int(1e3))
    n = len(m0);
    mout = np.zeros((n, niter));
    sigma = .5
    mout[:, 0] = m0[:,0]
    current = m0
    lMAP = -1*np.inf
    mMAP = current
    nacc = 0
    step = kwargs.get("step",.001)
    k = 0


    while k < niter:
        
        # generate a random guess (more/less...)
        candidate = Generate(current, step)
        
        # log-liklihood
        llcandidate = logliklihood(G, candidate, theta, sigma)

        # ------ For the Log Prior ----------#    
        uniform_current = [
                uniform(m0[0],0, 1),
                uniform(m0[1],0, 0.5),
                uniform(m0[2],0, 0.25),
                uniform(m0[3],0, 2.0),
                uniform(m0[4],0,.5)]

        
        uniform_candidate = [
                uniform(candidate[0],vmin,vmax),
                uniform(candidate[1],vmin,vmax),
                uniform(candidate[2],vmin,vmax),
                uniform(candidate[3],vmin,vmax),
                uniform(candidate[4],vmin,vmax)]
        
        # log prior current 
        lpcurrent = logprior(uniform_current)
        
        #log-prior 
        lpcandidate = logprior(uniform_candidate)
        # ------------------------------------# 
        
        # whatever these are called
        lr1 = logproposal(candidate, current, step)
        lr2 = logproposal(current, candidate, step)
        
        # more stuff

        llcurrent = logliklihood(G, current, theta, sigma)
        logalpha  = lpcandidate + llcandidate + lr1 - lpcurrent - llcurrent - lr2

        # logalpha gt. 0 
        if logalpha > 0:
            logalpha = 0
        else:
            pass

        # log of random normal 
        #logt = np.log(np.random.rand(1))
        
        # log of random uniform
        logt = np.log(np.random.uniform(low=0, high=1, size=1))

        # accept/reject the step
        if logt < logalpha:
            current = candidate
            nacc = nacc + 1        #number of 'acceptances'
            
            # update mMAP, lMAP
            if (lpcandidate + llcandidate) > lMAP:
                lMAP = lpcandidate + llcandidate
                mMAP = candidate
            else:
                pass
            # done 
       # 
        else:
            # do nothing 
            pass 
        mout[:, k] = current[:,0]
        accrate = np.float(nacc) / np.float(niter)
        k = k+1
    #end whie
    print(nacc,k)
    return mout 



def lag_corr(x,p):
    size = x.shape[1]
    x = x[p, :] # range(0,size,10)]
    corr = []
    step = 100
    lagl = []
    for lag in xrange(1, size-step, step):
        a = np.roll(x, shift=lag) 
        corr.append(np.corrcoef(a,x)[0,1])
        lagl.append(lag)
    return corr

def Fig117():
    fig, axx = plt.subplots(5,1)
    good_vals = 500
    for i,ax in enumerate(axx):

        series = mout[i, -good_vals:] 
        mn = np.mean(series)
        ax.scatter(range(good_vals), series)
        ax.axhline(m_real[i], color='red')
        ax.set_ylabel(r'$m_{}$'.format(i))
        if i == 4:
            ax.set_xlabel('sample number')
    pass 

def BigPlot():
    fig,ax = plt.subplots(5,5)
    fig.set_size_inches(8,8)
    for i in range(5):
        for j in range(5):
            if i == j:
                ax[i,j].hist(mout[i, -500:])
                ax[i,j].set_xticks([])
                ax[i,j].set_yticks([])
            else:
                xax = mout[i, -500:].copy()
                yax = mout[j, -500:].copy()
                xbnd = np.sort(xax)[int(500*.05)], np.sort(xax)[int(500*.95)]
                ybnd= np.sort(yax)[int(500*.05)], np.sort(yax)[int(500*.95)]
                
                # plot the mean()
                xmn = xax.mean()
                ymn = yax.mean()
                ax[i,j].scatter(xax, yax, color='red', alpha=.25)
                ax[i,j].scatter(xmn, ymn, marker = '+', color='blue')
                ax[i,j].axvline(xbnd[0])
                ax[i,j].axvline(xbnd[1])
                ax[i,j].axhline(ybnd[0])
                ax[i,j].axhline(ybnd[1])
    # done      
                
            if i == 4:
                ax[i,j].set_xlabel(r'$m_{}$'.format(j))
            if j == 0:
                ax[i,j].set_ylabel(r'$m_{}$'.format(i))
            if j != 0:
                ax[i,j].set_yticks([])
            if i != 4:
                ax[i,j].set_xticks([])


    


if __name__=='__main__':
#    mout = MC(step=.1, niter=int(1e6))
    mout = np.load('mout.npy')
    
    fig,axx = plt.subplots(5,1)
    for i,ax in enumerate(axx):
        corr = lag_corr(mout[:,range(0,int(1e6),100)], i)
        xax = np.array(range(0,len(corr)))*100.
        ax.plot(xax, corr)
        ax.set_ylabel(r'A($m_{}$)'.format(i))
    axx[-1].set_xlabel('lag')
   
