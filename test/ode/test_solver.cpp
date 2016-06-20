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

#include "gtest/gtest.h"

#include <cinder/common/vector.hpp>
#include <cinder/integrator/euler.hpp>
#include <cinder/integrator/midpoint.hpp>
#include <cinder/integrator/runge_kutta.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/solver.hpp>

namespace cinder {

namespace {
/**
 * ODE with discontinuities which -- if integrated -- will produce a sawtooth
 * function.
 */
struct SawtoothODE {
	struct State : public VectorBase<State, Real, 2> {
		using VectorBase<State, Real, 2>::VectorBase;

		constexpr static State scale() {
			return State({1.0, 1.0});
		}
	};

	template <typename State2, typename System>
	static void init(Time, const State2 &, const System &)
	{
	}

	static State s0() { return State({0, 1.0}); }

	static Time next_discontinuity(Time t)
	{
		// Change the direction all 10 milliseconds
		const Real interval = 10e-3;
		return Time::sec(std::ceil((t + Time(1)).sec() / interval) * interval);
	}

	template <typename State2, typename System>
	static void handle_discontinuity(Time, State2 &s, System &)
	{
		s[1] *= -1.0f;
	}

	template <typename State2, typename System>
	static void update(Time, State2 &, System &)
	{
		// Do nothing here
	}

	template <typename State2, typename System>
	static State df(const State2 &s, const System &)
	{
		// Linearly move into the direction indicated by the second state
		// component
		return State({s[1], 0.0});
	}
};

struct SawtoothTestRecorder {
	template <typename State, typename System>
	static void record(Time t, const State &s, const System &, bool)
	{
		const Real interval = 10e-3;
		const int i = std::ceil((t + Time(1)).sec() / interval);
		const float phase = i * interval - t.sec();
		EXPECT_NEAR((i % 2 == 1) ? interval - phase : phase, s[0], 1e-6);
	}
};

template <typename Integrator>
static void test_linear()
{
	SawtoothODE ode;
	Integrator integrator;
	SawtoothTestRecorder recorder;
	NullController controller;

	make_solver(ode, integrator, recorder, controller).solve(1_s, 1_ms);
}
}

TEST(solver, sawtooth_euler) { test_linear<EulerIntegrator>(); }

TEST(solver, sawtooth_midpoint) { test_linear<MidpointIntegrator>(); }

TEST(solver, sawtooth_runge_kutta) { test_linear<RungeKuttaIntegrator>(); }

TEST(solver, sawtooth_dormand_prince)
{
	test_linear<DormandPrinceIntegrator>();
}
}
