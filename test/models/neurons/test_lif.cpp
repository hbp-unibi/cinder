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

#include <cinder/models/synapses/cur_exp.hpp>
#include <cinder/models/synapses/cond_exp.hpp>
#include <cinder/models/neurons/lif.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/solver.hpp>
#include <cinder/ode/recorder.hpp>

namespace cinder {
TEST(lif, cur_exp)
{
	// Spike times computed using exact NEST simulation
	const std::vector<Time> expected_spikes{
	    34.99319742_ms, 46.87646214_ms,  56.90537176_ms,  205.43159846_ms,
	    208.4919938_ms, 211.37284922_ms, 215.56020005_ms, 224.43595096_ms};

	DormandPrinceIntegrator integrator;
	NullRecorder recorder;
	NullController controller;
	std::vector<Time> spikes;

	// Create a linear integrate and fire neuron with a current based synapse
	auto current_source = make_current_source(
	    CurExp(2_nA, 10_ms, {10_ms, 30_ms, 40_ms, 50_ms, 100_ms, 200_ms, 202_ms,
	                         204_ms, 206_ms, 208_ms}));
	auto neuron = make_neuron<LIF>(current_source,
	                               [&spikes](Time t) { spikes.push_back(t); });

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s);

	// Make sure the two results are within the range of one millisecond
	ASSERT_EQ(expected_spikes.size(), spikes.size());
	for (size_t i = 0; i < spikes.size(); i++) {
		EXPECT_GT(0.1_ms, std::abs(expected_spikes[i] - spikes[i]));
	}
}

TEST(lif, cond_exp)
{
	// Spike times computed using NEST simulation with 0.01 ms timestep
	const std::vector<Time> expected_spikes{
	    13.615_ms,  19.576_ms,  31.02_ms,   34.379_ms,  39.841_ms,  42.295_ms,
	    45.275_ms,  49.717_ms,  52.051_ms,  54.704_ms,  58.416_ms,  65.096_ms,
	    102.324_ms, 107.066_ms, 202.705_ms, 204.405_ms, 205.697_ms, 206.83_ms,
	    207.935_ms, 208.873_ms, 209.843_ms, 210.909_ms, 212.093_ms, 213.427_ms,
	    214.959_ms, 216.764_ms, 218.97_ms,  221.832_ms, 225.998_ms, 235.149_ms};

	DormandPrinceIntegrator integrator;
	NullRecorder recorder;
	NullController controller;
	std::vector<Time> spikes;

	// Create a linear integrate and fire neuron with a current based synapse
	auto current_source = make_current_source(
	    CondExp(0.1_uS, 10_ms, 0_V, {10_ms, 30_ms, 40_ms, 50_ms, 100_ms, 200_ms,
	                                 202_ms, 204_ms, 206_ms, 208_ms}));
	auto neuron = make_neuron<LIF>(current_source, [&spikes](Time t) {
		spikes.push_back(t);
	});

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s);

	// Make sure the two results are within the range of one millisecond
	ASSERT_EQ(expected_spikes.size(), spikes.size());
	for (size_t i = 0; i < spikes.size(); i++) {
		EXPECT_GT(0.3_ms, std::abs(expected_spikes[i] - spikes[i]));
	}
}
}
