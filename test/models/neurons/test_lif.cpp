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

#include <cinder/models/synapses/cur_alpha.hpp>
#include <cinder/models/synapses/cur_exp.hpp>
#include <cinder/models/synapses/cond_alpha.hpp>
#include <cinder/models/synapses/cond_exp.hpp>
#include <cinder/models/neurons/lif.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/solver.hpp>
#include <cinder/ode/recorder.hpp>

namespace cinder {
TEST(lif, cur_alpha)
{
	// Spike times computed using NEST simulation with 0.01ms timestep
	const std::vector<Time> expected_spikes{
	    22.57_ms,   33.225_ms,  39.564_ms,  44.511_ms,  48.701_ms,  52.747_ms,
	    56.239_ms,  59.812_ms,  63.856_ms,  68.956_ms,  76.934_ms,  106.55_ms,
	    115.764_ms, 207.214_ms, 209.748_ms, 211.704_ms, 213.487_ms, 215.216_ms,
	    216.952_ms, 218.739_ms, 220.618_ms, 222.633_ms, 224.843_ms, 227.332_ms,
	    230.236_ms, 233.805_ms, 238.618_ms, 246.892_ms};

	DormandPrinceIntegrator integrator;
	NullRecorder recorder;
	NullController controller;
	std::vector<Time> spikes;

	// Create a linear integrate and fire neuron with a current based synapse
	auto current_source = make_current_source(
	    CurAlpha(2_nA, 10_ms, {10_ms, 30_ms, 40_ms, 50_ms, 100_ms, 200_ms,
	                           202_ms, 204_ms, 206_ms, 208_ms}));
	auto neuron = make_neuron<LIF>(current_source,
	                               [&spikes](Time t) { spikes.push_back(t); });

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s);

	// Make sure the two results are within the range of one millisecond
	ASSERT_EQ(expected_spikes.size(), spikes.size());
	for (size_t i = 0; i < spikes.size(); i++) {
		EXPECT_GT(0.2_ms, std::abs(expected_spikes[i] - spikes[i]));
	}
}

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

TEST(lif, cond_alpha)
{
	// Spike times computed using NEST simulation with 0.01 ms timestep
	const std::vector<Time> expected_spikes{
	    15.756_ms,  18.864_ms,  21.808_ms,  24.883_ms,  28.321_ms,  31.85_ms,
	    34.248_ms,  36.347_ms,  38.37_ms,   40.411_ms,  42.181_ms,  43.736_ms,
	    45.204_ms,  46.64_ms,   48.079_ms,  49.545_ms,  51.013_ms,  52.346_ms,
	    53.605_ms,  54.829_ms,  56.043_ms,  57.265_ms,  58.509_ms,  59.788_ms,
	    61.116_ms,  62.508_ms,  63.982_ms,  65.56_ms,   67.271_ms,  69.155_ms,
	    71.271_ms,  73.71_ms,   76.628_ms,  80.338_ms,  85.651_ms,  99.951_ms,
	    104.989_ms, 107.927_ms, 110.676_ms, 113.5_ms,   116.58_ms,  120.145_ms,
	    124.62_ms,  131.219_ms, 204.507_ms, 206.17_ms,  207.385_ms, 208.397_ms,
	    209.274_ms, 210.068_ms, 210.812_ms, 211.523_ms, 212.212_ms, 212.886_ms,
	    213.55_ms,  214.208_ms, 214.863_ms, 215.518_ms, 216.175_ms, 216.836_ms,
	    217.503_ms, 218.178_ms, 218.863_ms, 219.559_ms, 220.269_ms, 220.995_ms,
	    221.738_ms, 222.502_ms, 223.289_ms, 224.102_ms, 224.944_ms, 225.819_ms,
	    226.732_ms, 227.688_ms, 228.694_ms, 229.757_ms, 230.887_ms, 232.096_ms,
	    233.4_ms,   234.82_ms,  236.383_ms, 238.129_ms, 240.117_ms, 242.441_ms,
	    245.265_ms, 248.922_ms, 254.324_ms};

	DormandPrinceIntegrator integrator;
	NullRecorder recorder;
	NullController controller;
	std::vector<Time> spikes;

	// Create a linear integrate and fire neuron with a current based synapse
	auto current_source = make_current_source(CondAlpha(
	    0.1_uS, 10_ms, 0_V, {10_ms, 30_ms, 40_ms, 50_ms, 100_ms, 200_ms, 202_ms,
	                         204_ms, 206_ms, 208_ms}));
	auto neuron = make_neuron<LIF>(current_source,
	                               [&spikes](Time t) { spikes.push_back(t); });

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s);

	// Make sure the two results are within the range of one millisecond
	ASSERT_EQ(expected_spikes.size(), spikes.size());
	for (size_t i = 0; i < spikes.size(); i++) {
		EXPECT_GT(0.5_ms, std::abs(expected_spikes[i] - spikes[i]));
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
	auto neuron = make_neuron<LIF>(current_source,
	                               [&spikes](Time t) { spikes.push_back(t); });

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s);

	// Make sure the two results are within the range of one millisecond
	ASSERT_EQ(expected_spikes.size(), spikes.size());
	for (size_t i = 0; i < spikes.size(); i++) {
		EXPECT_GT(0.21_ms, std::abs(expected_spikes[i] - spikes[i]));
	}
}
}
