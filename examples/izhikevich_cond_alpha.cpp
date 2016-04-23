/*
 *  Cinder -- C++ Single Spiking Neuron Simulator
 *  Copyright (C) 2015, 2016  Andreas Stöckel, Christoph Jenzen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file izhikevich_cond_alpha.cpp
 *
 * This example shows how to simulate a single neuron with two conductance based
 * synapses with alpha-function shaped conductivity trace. One of the synapses
 * is setup as an excitatory synapse, the other synapse is setup as an
 * inhibitory synapse.
 *
 * Note that the currents are injected into the Izhikevich neuron under the
 * assumption of a 1nF membrane capacitance. In the original Izhikevich model
 * the injected current has no particular unit -- the capacitance is needed to
 * convert between currents in ampere and a voltage differential in volts per
 * second.
 *
 * @author Andreas Stöckel
 */

#include <iostream>

#include <cinder/cinder.hpp>

using namespace cinder;

int main()
{
	// Use the DormandPrince integrator with default target error
	DormandPrinceIntegrator integrator;

	// Record the neuron state as CSV to cout
	CSVRecorder recorder(std::cout);

	// Use the NeuronController class to automatically abort the simulation
	// once the neuron has settled to its resting state
	NeuronController controller;

	// Assemble current source (here: two synapses)
	auto current_source = make_current_source(
	    CondAlpha(20_nS, 10_ms, 30_mV, {95_ms, 100_ms, 102_ms, 110_ms, 499_ms,
	                                     505_ms, 510_ms, 520_ms, 530_ms}),
	    CondAlpha(50_nS, 30_ms, -80_mV, {390_ms, 400_ms, 420_ms, 450_ms}));

	// Assemble the neuron (an Izhikevich neuron with default parameters)
	auto neuron = make_neuron<Izhikevich>(current_source, [] (Time) {
		// Lambda called whenever the neuron spikes
	});

	// Simulate the resulting ODE -- the recorder will record the results
	make_solver(neuron, integrator, recorder, controller).solve();
}

