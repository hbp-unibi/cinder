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

#include <cinder/models/synapses/cur_exp.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/controller.hpp>
#include <cinder/ode/solver.hpp>

namespace cinder {
namespace {

static const std::vector<Spike> input_spikes1 = {10_ms, 100_ms, 300_ms,
                                                 50_ms, 180_ms, 120_ms};
static const std::vector<Spike> input_spikes2 = {40_ms, 180_ms, 500_ms, 200_ms};
static const std::vector<Spike> input_spikes3 = {20_ms, 55_ms, 600_ms, 700_ms};

struct CurExpTestRecorder {
	Real calc_expected_current(Time t, Current cur, Time tau,
	                           const std::vector<Spike> &spikes)
	{
		Real res = 0.0;
		for (const Spike &spike : spikes) {
			Time t2 = t - spike.t;
			if (t2 >= 0_s) {
				res += cur.v() * std::exp(-t2.sec() / tau.sec());
			}
		}
		return res;
	}

	bool check_near_discontinuity(Time t, const std::vector<Spike> &spikes)
	{
		for (const Spike &spike : spikes) {
			if ((t - spike.t).abs() < 1_ms) {
				return true;
			}
		}
		return false;
	}

	template <typename State, typename System>
	void record(Time t, const State &s, const System &sys)
	{
		Real expected_current = 0.0;
		expected_current +=
		    calc_expected_current(t, 10_nA, 10_ms, input_spikes1);
		expected_current +=
		    calc_expected_current(t, 30_nA, 10_ms, input_spikes2);
		expected_current +=
		    calc_expected_current(t, -20_nA, 50_ms, input_spikes3);

		bool near_discontinuity = check_near_discontinuity(t, input_spikes1) ||
		                          check_near_discontinuity(t, input_spikes2) ||
		                          check_near_discontinuity(t, input_spikes3);
		if (!near_discontinuity) {
			EXPECT_NEAR(expected_current, sys.ode().current(s, sys), 1e-10);
		}
	}
};
}

TEST(cur_exp, basic)
{
	auto current_source =
	    make_current_source(CurExp(10_nA, 10_ms, input_spikes1),
	                        CurExp(30_nA, 10_ms, input_spikes2),
	                        CurExp(-20_nA, 50_ms));
	DormandPrinceIntegrator integrator;
	CurExpTestRecorder recorder;
	NullController controller;

	// Insert spikes into the exising synapse
	for (Spike spike: input_spikes3) {
		current_source.get<2>().input_spikes().push(spike);
	}

	make_solver(current_source, integrator, recorder, controller)
	    .solve(1_s, 1_ms);
}
}
