#!/usr/bin/env python

import sys
import pyNN.nest as sim
import time

sim.setup(timestep=0.001)

source = sim.Population(1, sim.SpikeSourceArray,
                        {"spike_times": [10, 11, 12, 13, 14, 15, 30, 32, 34,
                                         36, 38, 40]})

pop_izhikevich = sim.Population(1, sim.Izhikevich, {
    'a': 0.02,
    'b': 0.2,
    'c': -65.0,
    'd': 2.0
})
pop_izhikevich.record(["spikes"])

sim.Projection(source,
               pop_izhikevich,
               connector=sim.FromListConnector([(0, 0, 20.0, 0.1)]))

sim.run(1000.0)

print "Izhikevich neuron spike train"
print pop_izhikevich.get_data().segments[0].spiketrains

sim.end()

