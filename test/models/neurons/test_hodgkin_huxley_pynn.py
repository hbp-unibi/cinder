#!/usr/bin/env python

import sys
import numpy as np
import pyNN.nest as sim

timestep=0.1
sim.setup(timestep=timestep)

source = sim.Population(1, sim.SpikeSourceArray,
                        {"spike_times": [10, 30, 40, 50, 100, 200, 202, 204,
                                         206, 208]})

pop_cond_exp = sim.Population(1,
                              sim.HH_cond_exp,
                              {
                                  'gbar_Na': 20.0,
                                  'gbar_K': 6.0,
                                  'g_leak': 0.01,
                                  'cm': 0.2,
                                  'v_offset': -63.0,
                                  'e_rev_Na': 50.0,
                                  'e_rev_K': -90.0,
                                  'e_rev_leak': -65.0,
                                  'e_rev_E': 0.0,
                                  'e_rev_I': -80.0,
                                  'tau_syn_E': 10.0,
                                  'tau_syn_I': 2.0,
                                  'i_offset': 0.0,
                              })
pop_cond_exp.record(["spikes", "v"])

sim.Projection(source,
               pop_cond_exp,
               connector=sim.FromListConnector([(0, 0, 0.1, 0.1)]))

sim.run(1000.0)

print "Conductance based exponential HH neuron"
print pop_cond_exp.get_data().segments[0].spiketrains

voltage_trace = pop_cond_exp.get_data().segments[0].analogsignalarrays[0]
n = len(voltage_trace);
arr = np.zeros((n, 2));
arr[:, 0] = np.linspace(0, timestep * (n - 1), n)
arr[:, 1] = voltage_trace[:, 0]
np.savetxt("hh_trace_pynn.csv", arr, delimiter=",")


sim.end()

