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
#include <cinder/models/neurons/hodgkin_huxley.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/solver.hpp>
#include <cinder/ode/recorder.hpp>

namespace cinder {
TEST(hodgkin_huxley, cond_exp)
{
	// Spike times computed using NEST simulation with 0.1 ms timestep; Spike at
	// 209.2 is not emitted by the NEST simulation for unknown reasons
	const std::vector<Time> expected_spikes{
	    11.1_ms,  14.3_ms,  17.8_ms,  22._ms,   27.8_ms,  31.4_ms,  34.4_ms,
	    37.7_ms,  40.9_ms,  43.6_ms,  46.5_ms,  49.7_ms,  52.2_ms,  54.9_ms,
	    57.8_ms,  61.2_ms,  65.2_ms,  70.6_ms,  79.1_ms,  100.9_ms, 104._ms,
	    107.5_ms, 111.7_ms, 117.2_ms, 126.4_ms, 201.1_ms, 203.4_ms, 205.5_ms,
	    207.6_ms, 209.2_ms, 211.2_ms, 213.3_ms, 215.4_ms, 217.7_ms, 220.3_ms,
	    223.1_ms, 226.3_ms, 230.1_ms, 234.9_ms, 241.9_ms, 257.3_ms};

	DormandPrinceIntegrator integrator;
	NullRecorder recorder;
	NullController controller;
	std::vector<Time> spikes;

	// Create an AdEx neuron with a current based synapse
	auto current_source = make_current_source(
	    CondExp(0.1_uS, 10_ms, 0_V, {10_ms, 30_ms, 40_ms, 50_ms, 100_ms, 200_ms,
	                                 202_ms, 204_ms, 206_ms, 208_ms}));
	auto neuron = make_neuron<HodgkinHuxley>(
	    current_source, [&spikes](Time t) { spikes.push_back(t); });

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s);

	// Make sure the two results are within the range of one millisecond
	ASSERT_EQ(expected_spikes.size(), spikes.size());
	for (size_t i = 0; i < spikes.size(); i++) {
		EXPECT_GT(0.7_ms, std::abs(expected_spikes[i] - spikes[i]));
	}
}
}
