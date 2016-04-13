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

#include <cinder/models/current_source.hpp>
#include <cinder/integrator/euler.hpp>
#include <cinder/integrator/midpoint.hpp>
#include <cinder/integrator/runge_kutta.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/solver.hpp>

namespace cinder {
namespace {
struct NullCurrentTestRecorder {
	template <typename State, typename System>
	static void record(Time, const State &s, const System &sys)
	{
		EXPECT_EQ(0.0, sys.ode().current(s, sys));
	}
};

struct ConstantCurrentTestRecorder {
	template <typename State, typename System>
	static void record(Time, const State &s, const System &sys)
	{
		EXPECT_NEAR(10e-3, sys.ode().current(s, sys), 1e-6);
	}
};

struct StepCurrentTestRecorder {
	template <typename State, typename System>
	static void record(Time t, const State &s, const System &sys)
	{
		if (std::abs(t - 100_ms) > 0.1_ms && std::abs(t - 900_ms) > 0.1_ms) {
			EXPECT_NEAR((t > 100_ms && t < 900_ms) ? 10e-3 : 0.0,
			            sys.ode().current(s, sys), 1e-6);
		}
	}
};

struct MultiCurrentTestRecorder {
	template <typename State, typename System>
	static void record(Time t, const State &s, const System &sys)
	{
		if (std::abs(t - 100_ms) > 0.1_ms && std::abs(t - 900_ms) > 0.1_ms) {
			EXPECT_NEAR((t > 100_ms && t < 900_ms) ? 10e-3 : -10e-3,
			            sys.ode().current(s, sys), 1e-6);
		}
	}
};

template <typename Integrator>
static void test_null_current()
{
	NullCurrentSource ode;
	Integrator integrator;
	NullCurrentTestRecorder recorder;
	NullController controller;

	make_solver(ode, integrator, recorder, controller).solve(1_s, 1_ms);
}

template <typename Integrator>
static void test_constant_current()
{
	ConstantCurrentSource ode(10.0_mA);
	Integrator integrator;
	ConstantCurrentTestRecorder recorder;
	NullController controller;

	make_solver(ode, integrator, recorder, controller).solve(1_s, 1_ms);
}

template <typename Integrator>
static void test_step_current()
{
	StepCurrentSource ode(10.0_mA, 100_ms, 900_ms);
	Integrator integrator;
	StepCurrentTestRecorder recorder;
	NullController controller;

	make_solver(ode, integrator, recorder, controller).solve(1_s, 1_ms);
}

template <typename Integrator>
static void test_multi_current()
{
	ConstantCurrentSource constant(-10.0_mA);
	StepCurrentSource step(20.0_mA, 100_ms, 900_ms);
	auto ode = make_current_source(constant, step);
	Integrator integrator;
	MultiCurrentTestRecorder recorder;
	NullController controller;

	make_solver(ode, integrator, recorder, controller).solve(1_s, 1_ms);
}
}

TEST(current_source, null_current_euler)
{
	test_null_current<EulerIntegrator>();
}

TEST(current_source, null_current_midpoint)
{
	test_null_current<MidpointIntegrator>();
}

TEST(current_source, null_current_runge_kutta)
{
	test_null_current<RungeKuttaIntegrator>();
}

TEST(current_source, null_current_dormand_prince)
{
	test_null_current<DormandPrinceIntegrator>();
}

TEST(current_source, constant_current_euler)
{
	test_constant_current<EulerIntegrator>();
}

TEST(current_source, constant_current_midpoint)
{
	test_constant_current<MidpointIntegrator>();
}

TEST(current_source, constant_current_runge_kutta)
{
	test_constant_current<RungeKuttaIntegrator>();
}

TEST(current_source, constant_current_dormand_prince)
{
	test_constant_current<DormandPrinceIntegrator>();
}

TEST(current_source, step_current_euler)
{
	test_step_current<EulerIntegrator>();
}

TEST(current_source, step_current_midpoint)
{
	test_step_current<MidpointIntegrator>();
}

TEST(current_source, step_current_runge_kutta)
{
	test_step_current<RungeKuttaIntegrator>();
}

TEST(current_source, step_current_dormand_prince)
{
	test_step_current<DormandPrinceIntegrator>();
}

TEST(current_source, multi_current_euler)
{
	test_multi_current<EulerIntegrator>();
}

TEST(current_source, multi_current_midpoint)
{
	test_multi_current<MidpointIntegrator>();
}

TEST(current_source, multi_current_runge_kutta)
{
	test_multi_current<RungeKuttaIntegrator>();
}

TEST(current_source, multi_current_dormand_prince)
{
	test_multi_current<DormandPrinceIntegrator>();
}
}
