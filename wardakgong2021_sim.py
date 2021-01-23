"""Runs the network model on a single machine.

Parameter details can be found in the companion article. Note that the data
folders must be created before running this file.
"""

import sys

if len(sys.argv) < 3:
    print('Usage: python %s N_THREADS 100_ALPHA' % sys.argv[0])
    quit()


from wardakgong2021_networkmodel import network, networkclassical

keys = ['dt', 'delay', 'epsilon', 'tauMem', 'theta', 'Vr', 'timeout',
        'Jepsilon', 'datafolder', 'simtime', 'g', 'eta', 'orderCE', 'J',
        'nthreads', 'alpha']
def globals_kwargs():
    return {key: value for key, value in globals().items() if key in keys}

nthreads, alpha100 = map(int, sys.argv[1:3])
network_function = network if alpha100 != 200 else networkclassical
alpha = float(alpha100) / 100

dt = 0.1  # the resolution in ms
delay = 1.5  # synaptic delay in ms
epsilon = 0.1  # connection probability
tauMem = 20.0  # time constant of membrane potential in ms
theta = 20.0  # membrane threshold potential in mV
Vr = 10.0  # reset potential, "V_reset" defaults to 0.0

timeout = 60*60*1  # 60 * 60 * n_hours
Jepsilon = 1e-3  # correlation accuracy (to within epsilon)


## data12: brunel figure 7C (data11), 100s, alpha=1.2,1.99,2
## (g=5,nuext/nuthr=2) (CE=1000,CI=250,J=0.1mV,D=1.5ms)
#datafolder = 'data12'
#simtime = 10000.0 * .01  # Simulation time in ms
#g = 5.0  # ratio inhibitory weight/excitatory weight
#eta = 2.0  # external rate relative to threshold rate
#orderCE = 1000
#J = 0.1  # postsynaptic amplitude in mV

#network_function(**globals_kwargs())


# data13: biologically realistic (data9), 100s, alpha = 1.2,1.99,2
# (CE=4000,CI=1000,J=0.02,g=4,eta=nuext/nuthr=1)
datafolder = 'data13_lowstimulus3'
simtime = 10000.0 * 1  # Simulation time in ms
g = 4.0  # ratio inhibitory weight/excitatory weight
eta = 1.0 * 0.3  # external rate relative to threshold rate
orderCE = 4000
J = 0.02  # postsynaptic amplitude in mV

network_function(**globals_kwargs())
