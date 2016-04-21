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

#include <cmath>
#include <fstream>
#include <functional>

#include "gtest/gtest.h"

#include <cinder/models/synapses/delta.hpp>
#include <cinder/models/neurons/izhikevich.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/solver.hpp>
#include <cinder/ode/recorder.hpp>

namespace cinder {
TEST(izhikevich, delta)
{
	// Spike times computed using NEST simulation with 0.01ms timestep
	const std::vector<Time> expected_spikes{11.563_ms, 13.493_ms, 15.677_ms,
	                                        34.327_ms, 41.233_ms};

	DormandPrinceIntegrator integrator;
	std::ofstream os("izhikevich_delta.csv");
	CSVRecorder recorder(os);
	NullController controller;
	std::vector<Time> spikes;

	// Create a linear integrate and fire neuron with a current based synapse
	auto current_source = make_current_source(
	    Delta(20_mV, {10_ms, 11_ms, 12_ms, 13_ms, 14_ms, 15_ms, 30_ms, 32_ms,
	                  34_ms, 36_ms, 38_ms, 40_ms}));
	auto neuron = make_neuron<Izhikevich>(
	    current_source, [&spikes](Time t) { spikes.push_back(t); });

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s);

	// Make sure the two results are within the range of one millisecond
	ASSERT_EQ(expected_spikes.size(), spikes.size());
	for (size_t i = 0; i < spikes.size(); i++) {
		EXPECT_GT(0.21_ms, std::abs(expected_spikes[i] - spikes[i]));
	}
}
}
