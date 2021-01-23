"""Network output analysis.

Author: Asem Wardak

This file provides functions to analyse the simulation output and
generate some of the figures in the main article; the other figures are
produced in the Wolfram Language notebooks, found at `nb/*`.
"""

try:
    np
except NameError:
    import numpy as np
    import matplotlib.pyplot as plt
    plt.ioff()
    from scipy.signal import welch
    from scipy.optimize import curve_fit
    from scipy.stats import lognorm,variation
    
    from random import sample
    
    from gc import collect
    from functools import partial
    
    # define pool context manager to support the 'with' statement in Python 2
    # https://stackoverflow.com/a/54720031
    from sys import version_info
    from multiprocessing import Pool
    if version_info[0] == 2:
        from contextlib import contextmanager
        @contextmanager
        def Pool_(*args, **kwargs):
            pool = Pool(*args, **kwargs)
            yield pool
            pool.terminate()
    else:
        Pool_ = Pool
        
        

# Data loading
# spike file format: neuron_id, spike_time
# membrane potential file format: neuron_id, time, membrane_potential
        
def loadspikesworker(fname, Nmin, Nmax, i):
    s = [[] for _ in xrange(Nmax-Nmin)]
    with open(fname+"-%02d.gdf"%i,'r') as f:
        for l in f:
            l = l.strip().split('\t')
            N = int(l[0])-1
            if (N>=Nmin) and (N<Nmax):
                s[N-Nmin].append(float(l[1]))
    collect()
    return s
def loadspikes(fname, Nmin, Nmax):
    """Loads spikes from neurons Nmin (inclusive) to Nmax (exclusive).
    
    Returns: a list of neurons; each neuron being represented as a list of
             spike times
    """
    with Pool_(analysisnthreads) as p:
        return map(lambda *x: sum(x,[]),
                   *p.map(partial(loadspikesworker,fname,Nmin,Nmax),
                          xrange(nthreads)))
    
def loadspikesrasterworker(fname, Nmin, Nmax, i):
    s = []
    with open(fname+"-%02d.gdf"%i,'r') as f:
        for l in f:
            l = l.strip().split('\t')
            N = int(l[0])-1
            if (N>=Nmin) and (N<Nmax):
                s.append((N+1,float(l[1])))
    collect()
    return s
def loadspikesraster(fname, Nmin, Nmax):
    """Loads spikes from neurons Nmin (inclusive) to Nmax (exclusive).
    
    Returns: a list of pairs of neurons and spike times
    """
    with Pool_(analysisnthreads) as p:
        return np.array(sum(p.map(partial(loadspikesrasterworker,fname,Nmin,Nmax
                                          ),xrange(nthreads)),[]))

def loadspikesglobalworker(fname, Nmin, Nmax, i):
    s = []
    with open(fname+"-%02d.gdf"%i,'r') as f:
        for l in f:
            l = l.strip().split('\t')
            N = int(l[0])-1
            if (N>=Nmin) and (N<Nmax):
                s.append(float(l[1]))
    collect()
    return s
def loadspikesglobal(fname, Nmin, Nmax):
    """Loads a list of spike times occurring between neurons Nmin and Nmax."""
    with Pool_(analysisnthreads) as p:
        return sum(p.map(partial(loadspikesglobalworker,fname,Nmin,Nmax),
                         xrange(nthreads)), [])

def loadmultimetertimeworker(fname, n, i):
    s = np.zeros(T+1)
    with open(fname+"-%02d.dat"%i,'r') as f:
        for l in f:
            l = l.strip().split('\t')
            if int(float(l[0]))-1 == n:
                s[int(float(l[1]))] = float(l[2])
    collect()
    return s
def loadmultimetertime(fname, n):
    """Loads the membrane potential time series for neuron `n`."""
    with Pool_(analysisnthreads) as p:
        return sum(p.map(partial(loadmultimetertimeworker,fname,n),
                         xrange(nthreads)), np.zeros(T+1))

#def loadmultimetervaluesworker(fname, tmin, tmax, i):
#    s = []
#    with open(fname+"-%02d.dat"%i,'r') as f:
#        for l in f:
#            l = l.strip().split('\t')
#            t = int(float(l[1])  # in ms
#            if (t>=tmin) and (t<tmax):
#                s.append(float(l[2]))
#    collect()
#    return s
#def loadmultimetervalues(fname, tmin, tmax):
#   # load membrane potential values from time tmin to tmax
#    with Pool_(analysisnthreads) as p:
#        return sum(p.map(partial(loadmultimetervaluesworker,fname,tmin,tmax),
#                         xrange(nthreads)), [])
    
def loadglobalrateworker(fname, Ntot, i):
    """Loads a time series of the network firing rate."""
    s = np.zeros(T*10+1) # T/dt = 10s/0.1ms
    with open(fname+"-%02d.gdf"%i,'r') as f:
        for l in f:
            l = l.strip().split('\t')
            s[int(float(l[1])*10)] += 1 # times recorded in ms in file
    collect()
    return s*(10000./Ntot) # sampling rate: 10000 Hz
def loadglobalrate(fname, Ntot):
    with Pool_(analysisnthreads) as p:
        return sum(p.map(partial(loadglobalrateworker,fname,Ntot),
                         xrange(nthreads)),np.zeros(T*10+1))


# Plotting

def rasterplot(ax, spikesfilename, **kwargs):
    spikes = loadspikesraster(spikesfilename,0,50)
    ax.plot(spikes[:,1], spikes[:,0], '|', **kwargs)
    collect()
    
def Vtplot(ax, multimeterfilename, spikesfilename, mneuronid, sneuronid,
           **kwargs):
    Vt = loadmultimetertime(multimeterfilename,mneuronid)
    ax.plot(Vt,lw=1,**kwargs)
    for t in loadspikesglobal(spikesfilename,sneuronid,sneuronid+1):
        ax.vlines(t,Vt[int(t)],35,'r',lw=1)
    ax.set_ylim([-20,25])
    collect()
    
def nuFTplot(ax, spikesfilename, Ntot, **kwargs):
    globalrate = loadglobalrate(spikesfilename, Ntot)
    f,p = welch(globalrate,10000,nperseg=10000)
    ax.loglog(f,p,'.',**kwargs)
    collect()
    
def spectrogramplot(spikesfilename, Ntot, Tslice, **kwargs):
    globalrate = loadglobalrate(spikesfilename, Ntot)
    eng.cwt(matlab.double(globalrate[-Tslice*10:].tolist()),10000.,
            'VoicesPerOctave',48.,nargout=0)
    eng.colormap('jet')
    collect()
    
def nu0distplot(ax,spikesfilename,Nmin,Nmax):
    spikes = loadspikes(spikesfilename,Nmin,Nmax)
    bins = 10**np.linspace(-1,2.5,10)
    bincentres = (bins[:-1]+bins[1:])/2.
    nu0ss = np.transpose([np.histogram(x,bins=np.linspace(0,T,10))[0]/(T/1000.)
                          for x in spikes])
    hs = np.transpose([np.histogram(nu0s,bins=bins,density=True)[0]
                       for nu0s in nu0ss])
    # removing the first bin of half-line width
    bincentres = bincentres[1:]
    hs = hs[1:]
    ax.errorbar(bincentres,map(np.mean,hs),
                yerr=map(np.std,hs),fmt='.',capsize=2)
    # fitting
    xfits = 10**np.linspace(-1,2.5,100)
#    ax.plot(xfits,lognorm.pdf(xfits,*lognorm.fit(nu0s)))
    popt,_ = curve_fit(lambda x,s,scale: lognorm.pdf(x,s,0,scale),
                       bincentres,map(np.mean,hs))
    ax.plot(xfits,lognorm.pdf(xfits,popt[0],0,popt[1]),)
    ax.set_xscale("log")
    collect()
    
def CVISIplotworker(x):
    return variation(np.diff(x))
def CVISIplot(ax,spikesfilename,Nmin,Nmax):
    spikes = loadspikes(spikesfilename,Nmin,Nmax)
    with Pool_(20) as p:
        CVs = p.map(CVISIplotworker,spikes)
    CVs = [CV for CV in CVs if not np.isnan(CV)]
    ax.hist(CVs,bins=50,density=True,histtype='step')
    collect()

tau_on_dt = 20
def inputdistplot(ax1, ax2, multimeterfilename, mneuronid, **kwargs):
    binmin, binmax = -1, 3
    Vt = loadmultimetertime(multimeterfilename,mneuronid)
    It = tau_on_dt*np.diff(Vt) + Vt[:-1]
    It = It[Vt[1:]!=10]  # excludes resetting mechanism
    It_centred = It - np.mean(It)
#    ax.hist(It,bins=50,density=True,histtype='step')
    bins = np.flip(-10**np.linspace(binmin, binmax, 100))
    bincentres = (bins[:-1]+bins[1:])/2.
    hs = np.histogram(It_centred,bins=bins,density=True)[0]
#    print hs
    ax1.plot(bincentres[hs!=0],hs[hs!=0],'.')
    ax1.set_xscale("symlog")
    ax1.set_xlim([min(bincentres), max(bincentres)])
    ax1.set_yscale("log")
    bins = 10**np.linspace(binmin, binmax, 100)
    bincentres = (bins[:-1]+bins[1:])/2.
    hs = np.histogram(It_centred,bins=bins,density=True)[0]
    ax2.plot(bincentres[hs!=0],hs[hs!=0],'.')
    ax2.set_xscale("log")
    ax2.set_xlim([min(bincentres), max(bincentres)])
    ax2.set_yscale("log")
    collect()
    
def ffplot(ax, spikesfilename, Nmin, Nmax, **kwargs):
    spikes = loadspikes(spikesfilename,Nmin,Nmax)
    b = eng.logspace(*[matlab.double([i]) for i in [0,3.5,64]])[0]
    fanocalc = [eng.fanocalc(matlab.double(neuron),b)[0]
                for neuron in spikes if neuron]
    ax.loglog(b, map(np.nanmean,zip(*fanocalc)), **kwargs)
    collect()

def varmeancountplot(ax, spikesfilename, Nmin, Nmax):
    spikes = loadspikes(spikesfilename,Nmin,Nmax)
    spikecounts = sample([np.histogram(neuron,bins=T/1000,range=(0,T))[0]
                          for neuron in spikes],1000)
#    spikecounts = [np.histogram(neuron,bins=T/1000,range=(0,T))[0]
#                   for neuron in spikes]
    ax.loglog(map(np.mean,spikecounts),map(np.var,spikecounts),
              '.g',markersize=1)
    ax.loglog([1e-2,1e3],[1e-2,1e3],ls='--',c='.3')
    collect()


# Figures

def fig1(figname,xlim): # raster plot and membrane potential trace
    fig, ax = plt.subplots(nrows=2, ncols=len(alphas),
                           sharex='col', sharey='row',
                           figsize=plt.figaspect(1./len(alphas)), squeeze=False)
    fig.subplots_adjust(hspace=0.15,wspace=0)
    for i, alpha in enumerate(alphas):
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
#        ispikesfilename = "%s/alpha%.2fispikes-%d"%(datafolder,alpha,
#                                                    datanum(alpha)+2)
        rasterplot(ax[0,i], espikesfilename,
                   c='g' if alpha==alphas[-1] else 'b', ms=2)
        ax[0,i].set(title=r'$\alpha = %g$'%alpha, xlim=xlim)
        multimeterfilename = "%s/alpha%.2fV_m-%d"%(datafolder,alpha,
                                                   datanum(alpha)+3)
        Vtplot(ax[1,i], multimeterfilename, espikesfilename, 0, 0,
               c='g' if alpha==alphas[-1] else 'b')
        ax[1,i].set(xlabel='$t$ (ms)')
        ax[0,i].text(-0, 1.05, '(%s)'%'ab'[i], transform=ax[0,i].transAxes)
        #, size=20, weight='bold')
        ax[1,i].text(-0, 1.05, '(%s)'%'cd'[i], transform=ax[1,i].transAxes)
        #, size=20, weight='bold')
    for i, label in enumerate(['Neuron','$V$ (mV)']):
        ax[i,0].set(ylabel=label)
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)

def fig1_localisedinput(figname,xlim):
    alpha = alphas[0]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=plt.figaspect(0.125),
                           squeeze=False)
    espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,13502)
    multimeterfilename = "%s/alpha%.2fV_m-%d"%(datafolder,alpha,13503)
    Vtplot(ax[0,0], multimeterfilename, espikesfilename, 0, 0, c='purple')
    ax[0,0].set(title=r'$\alpha = %g$'%alpha, xlim=xlim)
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)
    
def specialfig1(figname,xlim,data_alphas):  
    fig, ax = plt.subplots(nrows=2, ncols=len(data_alphas),
                           sharex='col', sharey='row',
                           figsize=plt.figaspect(1./len(data_alphas)),
                           squeeze=False)
    fig.subplots_adjust(hspace=0.15,wspace=0)
    for i, j in enumerate(data_alphas):
        global datafolder, Nei_rec, nthreads, T, alphas
        datafolder, alpha = j
        Nei_rec, nthreads, T, alphas = data[datafolder]
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
        rasterplot(ax[0,i], espikesfilename,
                   c='g' if j==data_alphas[-1] else 'b', ms=2)
        ax[0,i].set(title=r'$\alpha = %g$'%alpha, xlim=xlim)
        multimeterfilename = "%s/alpha%.2fV_m-%d"%(datafolder,alpha,
                                                   datanum(alpha)+3)
        Vtplot(ax[1,i], multimeterfilename, espikesfilename, 0, 0,
               c='g' if j==data_alphas[-1] else 'b')
        ax[1,i].set(xlabel='$t$ (ms)')
        ax[0,i].text(-0, 1.05, '(%s)'%'ab'[i], transform=ax[0,i].transAxes)
        #, size=20, weight='bold')
        ax[1,i].text(-0, 1.05, '(%s)'%'cd'[i], transform=ax[1,i].transAxes)
        #, size=20, weight='bold')
    for i, label in enumerate(['Neuron','$V$ (mV)']):
        ax[i,0].set(ylabel=label)
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)

    
def fig2(figname): # power spectrum
    fig, ax = plt.subplots(nrows=1, ncols=len(alphas),
                           sharex='col', sharey='row',
                           figsize=plt.figaspect(1), squeeze=False)
#    fig.subplots_adjust(hspace=0,wspace=0)
    for i, alpha in enumerate(alphas):
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
#        ispikesfilename = "%s/alpha%.2fispikes-%d"%(datafolder,alpha,
#                                                    datanum(alpha)+2)
        nuFTplot(ax[0,i],espikesfilename,Nei_rec[0],
                 c='g' if alpha==alphas[-1] else 'b', ms=2)
#        ax[0,i].set(title=r'$\alpha = %g$'%alpha, xlim=[0.9,200],
#                    ylim=[1e-3,1e1], xlabel=r'$\omega_\nu$ (Hz)',
#                    ylabel='Intensity' if not i else '')
        ax[0,i].set(xlim=[0.9,200], ylim=[1e-3,1e1],
                    xlabel=r'$\omega_\nu$ (Hz)',
                    ylabel='Intensity' if not i else '')
    ax[0,0].loglog([10**1,10**2],[10**0,10**-1.5],c='tab:orange')
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)
    
def specialfig2(figname,data_alphas):
    fig, ax = plt.subplots(nrows=1, ncols=len(data_alphas),
                           sharex='col', sharey='row',
                           figsize=plt.figaspect(1), squeeze=False)
    for i, j in enumerate(data_alphas):
        global datafolder, Nei_rec, nthreads, T, alphas
        datafolder, alpha = j
        Nei_rec, nthreads, T, alphas = data[datafolder]
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
        nuFTplot(ax[0,i],espikesfilename,Nei_rec[0],
                 c='g' if j==data_alphas[-1] else 'b', ms=2)
        ax[0,i].set(xlim=[0.9,200], ylim=[1e-3,1e1],
                    xlabel=r'$\omega_\nu$ (Hz)',
                    ylabel='Intensity' if not i else '')
    ax[0,0].loglog([10**1,10**2],[10**0,10**-2],c='tab:orange')
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)
    
    
def fig3(figname,Tslice): # wavelet scalogram
    fig = eng.figure('Visible','off')
    for i, alpha in enumerate(alphas):
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
#        ispikesfilename = "%s/alpha%.2fispikes-%d"%(datafolder,alpha,
#                                                    datanum(alpha)+2)
        eng.subplot(1.,len(alphas)*1.,i+1.)
        spectrogramplot(espikesfilename,Nei_rec[0],Tslice)
        eng.title(r'\alpha=%g'%alpha)
        eng.ylim(matlab.double([np.log2(1e-3),np.log2(3e-1)]),nargout=0)
        eng.yticklabels(['%.1f'%(float(s)*1e3) if s else s
                         for s in eng.yticklabels()],nargout=0)
        eng.ylabel('Frequency (Hz)')
        eng.caxis(matlab.double([0.,[2.,20.][i]]),nargout=0)
        eng.ylabel(eng.colorbar(),'')
        eng.set(eng.gca(),'Position',
                eng.plus(eng.get(eng.gca(),'Position'),
                         matlab.double([-0.0,0.25,0.0,-0.25])))
        eng.set(eng.gca(),'TitleFontSizeMultiplier',2)
#        if i != len(alphas)-1: eng.colorbar('off',nargout=0)
        if i:
            eng.ylabel('',nargout=0)
            eng.yticklabels('',nargout=0)
            eng.set(eng.gca(),'Position',  # [-0.125,0,0.125,0]
                    eng.plus(eng.get(eng.gca(),'Position'),
                             matlab.double([-0.125,0.0,0.025,-0.0])))
#        eng.daspect(matlab.double([1,1,1]),nargout=0)
    eng.saveas(fig,figname,nargout=0)
    eng.close(fig,nargout=0)
    
def specialfig3(figname,Tslice,data_alphas):
    fig = eng.figure('Visible','off')
    for i, j in enumerate(data_alphas):
        global datafolder, Nei_rec, nthreads, T, alphas
        datafolder, alpha = j
        Nei_rec, nthreads, T, alphas = data[datafolder]
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
        eng.subplot(1.,len(data_alphas)*1.,i+1.)
        spectrogramplot(espikesfilename,Nei_rec[0],Tslice)
        eng.title(r'\alpha=%g'%alpha)
        eng.ylim(matlab.double([np.log2(1e-3),np.log2(3e-1)]),nargout=0)
        eng.yticklabels(['%.1f'%(float(s)*1e3) if s else s
                         for s in eng.yticklabels()],nargout=0)
        eng.ylabel('Frequency (Hz)')
        eng.caxis(matlab.double([0.,[2.,20.][i]]),nargout=0)
        eng.ylabel(eng.colorbar(),'')
        eng.set(eng.gca(),'Position',
                eng.plus(eng.get(eng.gca(),'Position'),
                         matlab.double([-0.0,0.25,0.0,-0.25])))
        eng.set(eng.gca(),'TitleFontSizeMultiplier',2)
        if i:
            eng.ylabel('',nargout=0)
            eng.yticklabels('',nargout=0)
            eng.set(eng.gca(),'Position',
                    eng.plus(eng.get(eng.gca(),'Position'),
                             matlab.double([-0.125,0.0,0.025,-0.0])))
    eng.saveas(fig,figname,nargout=0)
    eng.close(fig,nargout=0)
    
    
def fig4(figname): # network firing rate distribution
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex='col', sharey='row',
                           figsize=plt.figaspect(1), squeeze=False)
    alpha = alphas[0]
    espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                datanum(alpha)+1)
#    ispikesfilename = "%s/alpha%.2fispikes-%d"%(datafolder,alpha,
#                                               datanum(alpha)+2)
    #...(ispikesfilename,Nei_rec[0],sum(Nei_rec))
    nu0distplot(ax[0,0],espikesfilename,0,Nei_rec[0]) 
    ax[0,0].set(xlabel=r'Single-neuron $\nu_0$ (Hz)',ylabel='Density',
                title=r'$\alpha = %g$'%alpha)
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)
    
def fig5(figname): # CV of the ISI
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex='col', sharey='row',
                           figsize=plt.figaspect(1), squeeze=False)
    for alpha in [alphas[0]]:
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
    #    ispikesfilename = "%s/alpha%.2fispikes-%d"%(datafolder,alpha,
    #                                                datanum(alpha)+2)
        #...(ispikesfilename,Nei_rec[0],sum(Nei_rec))
        CVISIplot(ax[0,0],espikesfilename,0,Nei_rec[0]) 
        ax[0,0].set(xlabel=r'CV of ISI',ylabel='Density',
                    title=r'$\alpha = %g$'%alpha)
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)

def fig6(figname): # total input distribution to a given neuron
    nNeurons = 4
    fig, ax = plt.subplots(nrows=nNeurons, ncols=len(alphas)*2,
                           sharex='col', sharey='row',
                           figsize=plt.figaspect(0.5), squeeze=False)
    fig.subplots_adjust(hspace=0,wspace=0)
    for row in range(nNeurons):
        for i, alpha in enumerate(alphas):
            multimeterfilename = "%s/alpha%.2fV_m-%d"%(datafolder,alpha,
                                                       datanum(alpha)+3)
            inputdistplot(ax[row,2*i], ax[row,2*i+1], multimeterfilename, row)
            if not row: ax[row,2*i].set(title=r'$\alpha = %g$'%alpha)
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)

def fig7(figname):  # fano factor plot
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex='col', sharey='row',
                           figsize=plt.figaspect(1), squeeze=False)
    for alpha in [1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0]:
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
        ffplot(ax[0,0], espikesfilename, 0, Nei_rec[0],
               label=r'$\alpha = %g$'%alpha)
    ax[0,0].set(xlabel='T (ms)', ylabel='F(T)')
    ax[0,0].legend()
    ax[0,0].loglog([10**1,10**3],[10**0,10**1],ls='--',c='black')
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)

def specialfig7(figname,data_alphas):  # fano factor plot
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex='col', sharey='row',
                           figsize=plt.figaspect(1), squeeze=False)
    for i, j in enumerate(data_alphas):
        global datafolder, Nei_rec, nthreads, T, alphas
        datafolder, alpha = j
        Nei_rec, nthreads, T, alphas = data[datafolder]
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
        ffplot(ax[0,0], espikesfilename, 0, Nei_rec[0],
               label=r'$\alpha = %g$'%alpha)
    ax[0,0].set(xlabel='T (ms)', ylabel='F(T)')
    ax[0,0].legend()
    ax[0,0].loglog([10**1,10**3],[10**0,10**1],ls='--',c='black')
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)


def fig8(figname):  # spike count variance vs mean plot
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex='col', sharey='row',
                           figsize=plt.figaspect(1), squeeze=False)
    for alpha in [1.2]:
        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
                                                    datanum(alpha)+1)
        varmeancountplot(ax[0,0], espikesfilename, 0, Nei_rec[0])
        ax[0,0].set(title=r'$\alpha = %g$'%alpha,
                    xlim=[1e-2,1e3], ylim=[1e-2,1e3], xlabel='Mean (spikes)',
                    ylabel=r'Variance (spikes$^2$)')
    fig.savefig('%s'%figname, dpi=fig.dpi, bbox_inches='tight')
    plt.close(fig)
    


# Main
    
def datanum(alpha):
    if alpha == 2.00: return sum(Nei_rec)+1
    else: return sum(Nei_rec)+Nei_rec[0]/10
#    else: return sum(Nei_rec)+Nei_rec[0]/10 + 1


analysisnthreads = 6  # number of threads to perform data loading on
data = {'data11':([10000,2500],32,10000,[1.2,1.99]), #figure 7C (AI)
        'data12':([10000,2500],32,100000,[1.2,2.00]), #figure 7C long
        'data12_constantinput2':([10000,2500],32,10000,[1.2]),
#        'data12_localisedinput':([1,0],32,10000,[1.2]),
        'data12_reversalpotential':([10000,2500],32,10000,[1.2]),
        'data9':([40000,10000],32,10000,[1.2,1.99]), #biol. realistic (SI)
        'data13':([40000,10000],32,100000,[1.2,2.00]), #biol. realistic long
        'data13_constantinput2':([40000,10000],32,10000,[1.2]),
        }
        
        

# SI for alpha=1.2 and AI set for alpha=2
        
#collect()
#specialfig1('plots/specialfig1.pdf',[8050,8650],
#            [('data13',1.2),('data12',2.0)])

#collect()
#specialfig2('plots/specialfig2.pdf',[('data13', 1.2), ('data12', 2.0)])

#try:
#    eng
#except NameError:
#    import matlab.engine
#    eng = matlab.engine.start_matlab()
#collect()
#specialfig3('plots/specialfig3.pdf',10000,[('data13', 1.2), ('data12', 2.0)])

#try:
#    eng
#except NameError:
#    import matlab.engine
#    eng = matlab.engine.start_matlab()
#collect()
#specialfig7('plots/specialfig7.pdf',[('data13', 1.1), ('data13', 1.2),
#                                     ('data13', 1.4), ('data12', 2.0)])
        
        
        
        
        

for datafolder in ['data9']:#['data13']:#['data12']:#data:#.keys()[:1]:
    print datafolder,
    Nei_rec, nthreads, T, alphas = data[datafolder]
    
#    collect()
#    fig1('plots/%sfig1.pdf'%datafolder,[0,10000])#[8050,8650])
    

#    collect()
#    fig1_localisedinput('plots/%sfig1.pdf'%datafolder,[0,10000])
    
#    collect()
#    fig2('plots/%sfig2_2.pdf'%datafolder)
    
#    try:
#        eng
#    except NameError:
#        import matlab.engine
#        eng = matlab.engine.start_matlab()
#    collect()
#    fig3('plots/%sfig3_2.pdf'%datafolder,10000)
    
#    collect()
#    fig4('plots/%sfig4.pdf'%datafolder)
    
#    collect()
#    fig5('plots/%sfig5.pdf'%datafolder)

#    collect()
#    fig6('plots/%sfig6_2.pdf'%datafolder)

#    try:
#        eng
#    except NameError:
#        import matlab.engine
#        eng = matlab.engine.start_matlab()
#    collect()
#    fig7('plots/%sfig7.pdf'%datafolder)

    collect()
    fig8('plots/%sfig8.pdf'%datafolder)
    
#    alphas = [1.05,1.1,1.2,1.4,1.6,1.8,2]
#    alphas = [1.2]
#    for alpha in alphas:
#        espikesfilename = "%s/alpha%.2fespikes-%d"%(datafolder,alpha,
#                                                    datanum(alpha)+1)
#        espikes = loadspikes(espikesfilename,0,Nei_rec[0])
##        len_espikes = map(len,espikes)
#        # T is in ms
#        len_espikes = [len(x)/((T/1000)-0.002*len(x)) for x in espikes]  
#        print(np.mean(len_espikes), np.median(len_espikes),
#              np.percentile(len_espikes, [5, 95]), np.std(len_espikes))
#        plt.ion()
#        plt.figure()
#        plt.hist(len_espikes,100)
#        plt.title(alpha)
#        collect()
    




