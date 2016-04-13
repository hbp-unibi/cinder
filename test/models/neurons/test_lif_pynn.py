#!/usr/bin/env python

import sys
import pyNN.nest as sim

sim.setup(timestep=0.001);

source = sim.Population(1, sim.SpikeSourceArray,
                        {"spike_times": [10, 30, 40, 50, 100, 200, 202, 204,
                                         206, 208]})

pop_curr_exp = sim.Population(1, sim.IF_curr_exp, {
    "cm": 1.0,
    "tau_m": 20.0,
    "v_thresh": -50.0,
    "v_rest": -65.0,
    "v_reset": -65.0,
    "tau_refrac": 0.1,
    "tau_syn_E": 10.0
})
pop_curr_exp.record(["spikes"])

pop_cond_exp = sim.Population(1, sim.IF_cond_exp, {
    "cm": 1.0,
    "tau_m": 20.0,
    "v_thresh": -50.0,
    "v_rest": -65.0,
    "v_reset": -65.0,
    "tau_refrac": 0.1,
    "e_rev_E": 0.0,
    "tau_syn_E": 10.0
})
pop_cond_exp.record(["spikes"])

sim.Projection(source,
               pop_curr_exp,
               connector=sim.FromListConnector([(0, 0, 2.0, 0.1)]))
sim.Projection(source,
               pop_cond_exp,
               connector=sim.FromListConnector([(0, 0, 0.1, 0.1)]))

sim.run(1000.0);

print "Current based LIF neuron"
print pop_curr_exp.get_data().segments[0].spiketrains

print "Conductance based LIF neuron"
print pop_cond_exp.get_data().segments[0].spiketrains

sim.end();

