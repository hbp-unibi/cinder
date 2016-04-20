#!/usr/bin/env python

import sys
import pyNN.nest as sim
import time

sim.setup(timestep=0.001)

source = sim.Population(1, sim.SpikeSourceArray,
                        {"spike_times": [10, 30, 40, 50, 100, 200, 202, 204,
                                         206, 208]})

pop_cond_exp = sim.Population(1, sim.EIF_cond_exp_isfa_ista,
                              {
    'cm': 0.281,
    'tau_m': 9.3667,
    'v_thresh': -50.4,
    'v_rest': -70.6,
    'v_reset': -70.6,
    'v_spike': -40.0,
    'tau_refrac': 0.1,
    'delta_T': 2.0,
    'a': 4.0,
    'b': 0.0805,
    'tau_w': 144.0,
    'e_rev_E': 0.0,
    'tau_syn_E': 10.0,
})
pop_cond_exp.record(["spikes"])

sim.Projection(source,
               pop_cond_exp,
               connector=sim.FromListConnector([(0, 0, 0.1, 0.1)]))

sim.run(1000.0)

print "Conductance based AdEx neuron"
print pop_cond_exp.get_data().segments[0].spiketrains

sim.end()
