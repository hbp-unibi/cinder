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
 * @file izhikevich_sweep_constant_current.cpp
 *
 * Example program which injects an increasing current into an Izhikevich neuron
 * with default parameters and counts the number of output spikes. Outputs the
 * output spike count over the current.
 *
 * @author Andreas Stöckel
 */

#include <iostream>

#include <cinder/cinder.hpp>

using namespace cinder;

int main()
{
	// Write the CSV header
	std::cout << "current, spike_count" << std::endl;

	// Sweep over a small current range in a hundred steps
	for (Current i = 0_uA; i < 10_nA; i += 100_pA) {
		// Use the DormandPrince integrator with default target error
		DormandPrinceIntegrator integrator;

		// Do not record
		NullRecorder recorder;

		// Use the NeuronController class to automatically abort the simulation
		// once the neuron has settled to its resting state
		NullController controller;

		// Assemble the neuron (an Izhikevich neuron with default parameters)
		size_t spike_count = 0;
		auto neuron = make_neuron<Izhikevich>(
		    ConstantCurrentSource(i), [&spike_count](Time) { spike_count++; });

		// Simulate the resulting ODE -- the recorder will record the results
		make_solver(neuron, integrator, recorder, controller).solve(1_s);

		// Write the result of the exploration
		std::cout << i << ", " << spike_count << '\n';
	}
}

