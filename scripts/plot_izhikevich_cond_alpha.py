#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  Cinder -- C++ Single Spiking Neuron Simulator
#  Copyright (C) 2015, 2016  Andreas Stoeckel, Christoph Jenzen
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Script for plotting the results of the izhikevich_cond_alpha example program.
"""

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    print "Usage: " + sys.argv[0] + " <INPUT CSV>"
    sys.exit(1)


def cm2inch(value):
    return value / 2.54


data = np.loadtxt(sys.argv[1], delimiter=",")

fig = plt.figure(figsize=(cm2inch(21), cm2inch(11)))
ax = fig.add_subplot(3, 1, 1)
ax.plot(data[:, 0] * 1e3, data[:, 1] * 1e3, color='#000000')
ax.set_xlabel("Time [ms]")
ax.set_ylabel("Voltage [mV]")

ax = fig.add_subplot(3, 1, 2)
ax.plot(data[:, 0] * 1e3, data[:, 3] * 1e6, color='#000000')
ax.set_xlabel("Time [ms]")
ax.set_ylabel("Exc. Cond. [uS]")

ax = fig.add_subplot(3, 1, 3)
ax.plot(data[:, 0] * 1e3, data[:, 5] * 1e6, color='#000000')
ax.set_xlabel("Time [ms]")
ax.set_ylabel("Inh. Cond. [uS]")

fig.savefig("izhikevich_cond_alpha.pdf", format='pdf', bbox_inches='tight')

