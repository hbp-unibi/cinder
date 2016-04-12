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
	DormandPrinceIntegrator integrator;
	NullRecorder recorder;
	NullController controller;
	std::vector<Time> spikes;

	// Create a linear integrate and fire neuron with a current based synapse
	auto current_source = make_current_source(
	    CurExp(20_nA, 10_ms, {10_ms, 30_ms, 40_ms, 50_ms, 100_ms, 200_ms,
	                          202_ms, 204_ms, 206_ms, 208_ms}));
	auto neuron = make_neuron<LIF>(current_source,
	                               [&spikes](Time t) { spikes.push_back(t); });

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s, 1_ms);
}

TEST(lif, cond_exp)
{
	DormandPrinceIntegrator integrator;
	NullRecorder recorder;
	NullController controller;
	std::vector<Time> spikes;

	// Create a linear integrate and fire neuron with a current based synapse
	auto current_source = make_current_source(
	    CondExp(0.5_uS, 5_ms, 0_V, {10_ms, 30_ms, 40_ms, 50_ms, 100_ms, 200_ms,
	                          202_ms, 204_ms, 206_ms, 208_ms}));
	auto neuron = make_neuron<LIF>(current_source,
	                               [&spikes](Time t) { spikes.push_back(t); });

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s, 1_ms);
}


}
