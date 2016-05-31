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
#include <cinder/models/neurons/adex.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/solver.hpp>
#include <cinder/ode/recorder.hpp>

namespace cinder {
TEST(adex, cond_exp)
{
	// Spike times computed using NEST simulation with 0.01 ms timestep
	const std::vector<Time> expected_spikes{
	    11.894_ms,  14.12_ms,   17.115_ms,  22.074_ms,  30.969_ms,  32.832_ms,
	    35.218_ms,  38.681_ms,  41.068_ms,  42.615_ms,  44.504_ms,  46.987_ms,
	    50.435_ms,  51.823_ms,  53.482_ms,  55.582_ms,  58.579_ms,  103.466_ms,
	    107.495_ms, 202.538_ms, 203.64_ms,  204.595_ms, 205.419_ms, 206.283_ms,
	    206.953_ms, 207.67_ms,  208.371_ms, 208.974_ms, 209.614_ms, 210.296_ms,
	    211.028_ms, 211.819_ms, 212.682_ms, 213.634_ms, 214.699_ms, 215.916_ms,
	    217.348_ms, 219.115_ms, 221.513_ms, 226.102_ms};

	DormandPrinceIntegrator integrator;
	NullRecorder recorder;
	AutoController controller;
	std::vector<Time> spikes;

	// Create an AdEx neuron with a current based synapse
	auto current_source = make_current_source(
	    CondExp(0.1_uS, 10_ms, 0_V, {10_ms, 30_ms, 40_ms, 50_ms, 100_ms, 200_ms,
	                                 202_ms, 204_ms, 206_ms, 208_ms}));
	auto neuron = make_neuron<AdEx>(current_source, [&spikes](Time t) {
		spikes.push_back(t);
	}, AdEx::Parameters().v_reset(-70.6_mV).v_spike(-40.0_mV));

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve();

	// Make sure the two results are within the range of one millisecond
	ASSERT_EQ(expected_spikes.size(), spikes.size());
	for (size_t i = 0; i < spikes.size(); i++) {
		EXPECT_GT(0.3_ms, std::abs(expected_spikes[i] - spikes[i]));
	}
}
}
