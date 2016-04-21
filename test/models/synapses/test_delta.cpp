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
#include <functional>

#include "gtest/gtest.h"

#include <cinder/models/synapses/delta.hpp>
#include <cinder/models/neurons/lif.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/controller.hpp>
#include <cinder/ode/solver.hpp>

namespace cinder {
namespace {

static const std::vector<Spike> input_spikes = {10_ms, 100_ms, 300_ms};
using ODE = ODEBase<LIFState>;

struct DeltaTestRecorder {
	Real calc_expected_voltage(Time t, Voltage pulse_v,
	                           const std::vector<Spike> &spikes)
	{
		Real res = 0.0;
		for (const Spike &spike : spikes) {
			Time t2 = t - spike.t;
			if (t2 >= 0_s) {
				res += pulse_v * spike.w;
			}
		}
		return res;
	}

	bool check_near_discontinuity(Time t, const std::vector<Spike> &spikes)
	{
		for (const Spike &spike : spikes) {
			if (std::abs(t - spike.t) < 1_ms) {
				return true;
			}
		}
		return false;
	}

	template <typename State, typename System>
	void record(Time t, const State &s, const System &)
	{
		Real expected_current = 0.0;
		expected_current +=
		    calc_expected_voltage(t, 0.5_V, input_spikes);

		bool near_discontinuity = check_near_discontinuity(t, input_spikes);
		if (!near_discontinuity) {
			EXPECT_NEAR(expected_current, s[0], 1e-10);
		}
	}
};
}

TEST(delta, basic)
{
	MultiODE<ODE, Delta> ode({}, Delta(0.5_V, input_spikes));

	DormandPrinceIntegrator integrator;
	DeltaTestRecorder recorder;
	NullController controller;

	make_solver(ode, integrator, recorder, controller)
	    .solve(1_s, 1_ms);
}
}
