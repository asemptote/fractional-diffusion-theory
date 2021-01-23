"""Network model.

Author: Asem Wardak

This file provides functions containing simulation code for the classical,
homogeneous and fractional, heterogeneous network models.

Add the following to get

- constant input:
    ```
    constant_input = nest.Create("step_current_generator")
    nest.SetStatus(constant_input, {"amplitude_times": [dt, 5000.0],
                                    "amplitude_values": [0.0, CE*nu_ex * 10]})
    nest.Connect(constant_input, nodes_ex+nodes_in,
                 syn_spec={"model": "excitatory", "weight": J_ex})
    ```
  and modify the noise accordingly (or set
  `{"amplitude_times": [dt, 2*dt], "amplitude_values": [0.0, CE*nu_ex * 1]}`
  and remove the noise completely)

- reversal potential: add `"V_min": -10.0` to neuron_params

The variable parameters are defined in wardakgong2021_sim.py; these
functions are called from there.    
"""


import nest
import numpy
import sys
from time import time as tim
from numpy.random import pareto,seed
from multiprocessing import Pool, Value
from functools import partial
from scipy.special import gamma
# define pool context manager to support the 'with' statement in Python 2
# https://stackoverflow.com/a/54720031
if sys.version_info[0] == 2:
    from contextlib import contextmanager
    @contextmanager
    def Pool_(*args, **kwargs):
        pool = Pool(*args, **kwargs)
        yield pool
        pool.terminate()
else:
    Pool_ = Pool
    
nfailures_tightbalance = Value('i', 0)

# generate excitatory-weight-sum-correlated inhibitory weights for a neuron
def f(alpha,x0fac,x1fac,CE,CI,sum_a):
    seed()
    CEonCI = CE/CI
    epsilonXCE = Jepsilon*CE
    b = x1fac*(pareto(alpha,CI) + 1) + x0fac
    tf0 = tim()
    # mean of ensembles a,b is CI,CE
    while (abs(sum_a-CEonCI*sum(b)) > epsilonXCE) and (tim()-tf0 < timeout):
        b1 = x1fac*(pareto(alpha,CI) + 1) + x0fac
        if abs(sum_a-CEonCI*sum(b1)) < abs(sum_a-CEonCI*sum(b)):
            b = b1
#    print sum_a, sum(b), abs(sum_a-CEonCI*sum(b))/abs(sum_a)
    if tim() - tf0 > timeout:
        with nfailures_tightbalance.get_lock():
            nfailures_tightbalance.value += 1
    return b

def network(**kwargs):
    with nfailures_tightbalance.get_lock(): nfailures_tightbalance.value = 0
    globals().update(kwargs)

    nest.ResetKernel()
    startbuild = tim()
    order = int(orderCE / (epsilon * 4)) #2500
    NE = 4 * order  # number of excitatory neurons
    NI = 1 * order  # number of inhibitory neurons
    N_neurons = NE + NI  # number of neurons in total
#    N_rec = 50  # record from 50 neurons
    CE = int(epsilon * NE)  # number of excitatory synapses per neuron
    CI = int(epsilon * NI)  # number of inhibitory synapses per neuron
    #C_tot = int(CI + CE)  # total number of synapses per neuron
    neuron_params = {"C_m": 1.0, "tau_m": tauMem, "t_ref": 2.0, "E_L": 0.0,
                     "V_reset": Vr, "V_m": 0.0, "V_th": theta}
    J_ex = J  # amplitude of excitatory postsynaptic potential
    J_in = -g * J_ex  # amplitude of inhibitory postsynaptic potential
    nu_th = theta / (J * CE * tauMem)
    nu_ex = eta * nu_th
    p_rate = 1000.0 * nu_ex * CE
    nest.SetKernelStatus({"resolution": dt, "print_time": True,
                          "overwrite_files": True,
                          "local_num_threads": nthreads})
    nest.SetDefaults("iaf_psc_delta", neuron_params)
    nest.SetDefaults("poisson_generator", {"rate": p_rate/CE})
    nodes_ex = nest.Create("iaf_psc_delta", NE)
    nodes_in = nest.Create("iaf_psc_delta", NI)
    noise = nest.Create("poisson_generator", CE)
    espikes = nest.Create("spike_detector")
    ispikes = nest.Create("spike_detector")
    nest.SetStatus(espikes, [{"label": "%s/alpha%.2fespikes"%(datafolder,alpha),
                              "withtime": True, "withgid": True,
                              "to_file": True}])
    nest.SetStatus(ispikes, [{"label": "%s/alpha%.2fispikes"%(datafolder,alpha),
                              "withtime": True, "withgid": True,
                              "to_file": True}])
    nest.CopyModel("static_synapse", "excitatory",
                   {"weight": J_ex, "delay": delay})
    nest.CopyModel("static_synapse", "inhibitory",
                   {"weight": J_in, "delay": delay})
    A_alpha = gamma(1+alpha)*numpy.sin(numpy.pi*alpha/2)/numpy.pi
    D = 0.5
    # pareto pdf = alpha*x1**alpha/(x-x0)**(alpha+1), defined for x > x0+x1
    x1fac = (2*A_alpha*D/alpha)**(1/alpha)
    x0fac = 1 - x1fac*alpha/(alpha-1)
    J_noise_ex = J_ex * (x1fac*(pareto(alpha,(NE,CE)) + 1) + x0fac)
    J_noise_in = J_ex * (x1fac*(pareto(alpha,(NI,CE)) + 1) + x0fac)
    # correlated amplitude populations:
    samples_ex = x1fac*(pareto(alpha,(NE+NI,CE)) + 1) + x0fac
#    print x0fac,x1fac
   
    with Pool_(nthreads) as p:
        samples_in = numpy.array(p.map(partial(f,alpha,x0fac,x1fac,CE,CI),
                                       numpy.sum(samples_ex,1),1))
    J_ex_tot = J_ex * samples_ex
    J_in_tot = J_in * samples_in
    multimeter = nest.Create("multimeter")
    nest.SetStatus(multimeter, {"to_memory": False, "withtime": True,
                                "record_from": ["V_m"], "to_file": True,
                                "label": "%s/alpha%.2fV_m"%(datafolder,alpha)})
                                #"interval": 100.0,
    nest.Connect(multimeter,nodes_ex+nodes_in) # nodes_ex[:N_rec]+...
    nest.Connect(noise, nodes_ex,
                 syn_spec={"model": "excitatory", "weight": J_noise_ex})
    nest.Connect(noise, nodes_in,
                 syn_spec={"model": "excitatory", "weight": J_noise_in})
    nest.Connect(nodes_ex, espikes, syn_spec="excitatory") # nodes_ex[:N_rec]
    nest.Connect(nodes_in, ispikes, syn_spec="excitatory") # nodes_in[:N_rec]
    conn_params_ex = {'rule': 'fixed_indegree', 'indegree': CE}
    nest.Connect(nodes_ex, nodes_ex + nodes_in, conn_params_ex,
                 syn_spec={"model": "excitatory", "weight": J_ex_tot})
    conn_params_in = {'rule': 'fixed_indegree', 'indegree': CI}
    nest.Connect(nodes_in, nodes_ex + nodes_in, conn_params_in,
                 syn_spec={"model": "inhibitory", "weight": J_in_tot})
    endbuild = tim()
    nest.Simulate(simtime)
    endsimulate = tim()
    events_ex = nest.GetStatus(espikes, "n_events")[0]
    events_in = nest.GetStatus(ispikes, "n_events")[0]
    rate_ex = events_ex / simtime * 1000.0 / NE
    rate_in = events_in / simtime * 1000.0 / NI
    num_synapses = (nest.GetDefaults("excitatory")["num_connections"]
                    + nest.GetDefaults("inhibitory")["num_connections"])
    build_time = endbuild - startbuild
    sim_time = endsimulate - endbuild
    print("Number of tight balance failures: {0}"
          .format(nfailures_tightbalance.value))
    print("alpha             : {0}".format(alpha))
    print("Number of neurons : {0}".format(N_neurons))
    print("Number of synapses: {0}".format(num_synapses))
    print("       Exitatory  : {0}".format(int(CE * N_neurons) + N_neurons))
    print("       Inhibitory : {0}".format(int(CI * N_neurons)))
    print("Excitatory rate   : %.2f Hz" % rate_ex)
    print("Inhibitory rate   : %.2f Hz" % rate_in)
    print("Building time     : %.2f s" % build_time)
    print("Simulation time   : %.2f s" % sim_time)


def networkclassical(**kwargs):
    globals().update(kwargs)

    nest.ResetKernel()
    startbuild = tim()
    order = int(orderCE / (epsilon * 4)) #2500
    NE = 4 * order  # number of excitatory neurons
    NI = 1 * order  # number of inhibitory neurons
    N_neurons = NE + NI  # number of neurons in total
#    N_rec = 50  # record from 50 neurons
    CE = int(epsilon * NE)  # number of excitatory synapses per neuron
    CI = int(epsilon * NI)  # number of inhibitory synapses per neuron
    #C_tot = int(CI + CE)  # total number of synapses per neuron
    neuron_params = {"C_m": 1.0, "tau_m": tauMem, "t_ref": 2.0, "E_L": 0.0,
                     "V_reset": Vr, "V_m": 0.0, "V_th": theta}
    J_ex = J  # amplitude of excitatory postsynaptic potential
    J_in = -g * J_ex  # amplitude of inhibitory postsynaptic potential
    nu_th = theta / (J * CE * tauMem)
    nu_ex = eta * nu_th
    p_rate = 1000.0 * nu_ex * CE
    nest.SetKernelStatus({"resolution": dt, "print_time": True,
                          "overwrite_files": True,
                          "local_num_threads": nthreads})
    nest.SetDefaults("iaf_psc_delta", neuron_params)
    nest.SetDefaults("poisson_generator", {"rate": p_rate}) # `CE` gen's in 1
    nodes_ex = nest.Create("iaf_psc_delta", NE)
    nodes_in = nest.Create("iaf_psc_delta", NI)
    noise = nest.Create("poisson_generator")#, CE) in the superdiffusive case
    espikes = nest.Create("spike_detector")
    ispikes = nest.Create("spike_detector")
    nest.SetStatus(espikes, [{"label": "%s/alpha%.2fespikes"%(datafolder,alpha),
                              "withtime": True, "withgid": True,
                              "to_file": True}])
    nest.SetStatus(ispikes, [{"label": "%s/alpha%.2fispikes"%(datafolder,alpha),
                              "withtime": True, "withgid": True,
                              "to_file": True}])
    nest.CopyModel("static_synapse", "excitatory",
                   {"weight": J_ex, "delay": delay})
    nest.CopyModel("static_synapse", "inhibitory",
                   {"weight": J_in, "delay": delay})
    multimeter = nest.Create("multimeter")
    nest.SetStatus(multimeter, {"to_memory": False, "withtime": True,
                                "record_from": ["V_m"], "to_file": True,
                                "label": "%s/alpha%.2fV_m"%(datafolder,alpha)})
                                #"interval": 100.0,
    nest.Connect(multimeter,nodes_ex+nodes_in) # nodes_ex[:N_rec]+...
    nest.Connect(noise, nodes_ex, syn_spec="excitatory")
    nest.Connect(noise, nodes_in, syn_spec="excitatory")
    nest.Connect(nodes_ex, espikes, syn_spec="excitatory") # nodes_ex[:N_rec]
    nest.Connect(nodes_in, ispikes, syn_spec="excitatory") # nodes_in[:N_rec]
    conn_params_ex = {'rule': 'fixed_indegree', 'indegree': CE}
    nest.Connect(nodes_ex, nodes_ex + nodes_in, conn_params_ex, "excitatory")
    conn_params_in = {'rule': 'fixed_indegree', 'indegree': CI}
    nest.Connect(nodes_in, nodes_ex + nodes_in, conn_params_in, "inhibitory")
    endbuild = tim()
    nest.Simulate(simtime)
    endsimulate = tim()
    events_ex = nest.GetStatus(espikes, "n_events")[0]
    events_in = nest.GetStatus(ispikes, "n_events")[0]
    rate_ex = events_ex / simtime * 1000.0 / NE
    rate_in = events_in / simtime * 1000.0 / NI
    num_synapses = (nest.GetDefaults("excitatory")["num_connections"]
                    + nest.GetDefaults("inhibitory")["num_connections"])
    build_time = endbuild - startbuild
    sim_time = endsimulate - endbuild
    print("alpha             : {0}".format(alpha))
    print("Number of neurons : {0}".format(N_neurons))
    print("Number of synapses: {0}".format(num_synapses))
    print("       Exitatory  : {0}".format(int(CE * N_neurons) + N_neurons))
    print("       Inhibitory : {0}".format(int(CI * N_neurons)))
    print("Excitatory rate   : %.2f Hz" % rate_ex)
    print("Inhibitory rate   : %.2f Hz" % rate_in)
    print("Building time     : %.2f s" % build_time)
    print("Simulation time   : %.2f s" % sim_time)





