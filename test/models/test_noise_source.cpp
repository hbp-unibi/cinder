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

#include <cinder/models/noise_source.hpp>
#include <cinder/integrator/euler.hpp>
#include <cinder/integrator/midpoint.hpp>
#include <cinder/integrator/runge_kutta.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/recorder.hpp>
#include <cinder/ode/solver.hpp>

namespace cinder {
namespace {
struct GaussianNoiseSourceTestRecorder {
	std::vector<Real> is;

	template <typename State, typename System>
	void record(Time, const State &s, const System &sys, bool)
	{
		is.emplace_back(sys.ode().current(s, sys));
	}
};

template <typename Integrator>
static void test_gaussian_noise_current_source()
{
	GaussianNoiseSource ode(GaussianNoiseSourceParameters().i_stddev(1_A));
	Integrator integrator;
	GaussianNoiseSourceTestRecorder recorder;
	NullController controller;

	make_solver(ode, integrator, recorder, controller).solve(10_s, 1_ms);

	// Make sure the average is close to zero
	double avg = std::accumulate(recorder.is.begin(), recorder.is.end(), 0.0) /
	             double(recorder.is.size());
	EXPECT_NEAR(0.0, avg, 1e-4);
}
}

TEST(noise_source, gaussian_noise_current_source_euler)
{
	test_gaussian_noise_current_source<EulerIntegrator>();
}

TEST(noise_source, gaussian_noise_current_source_midpoint)
{
	test_gaussian_noise_current_source<MidpointIntegrator>();
}

TEST(noise_source, gaussian_noise_current_source_runge_kutta)
{
	test_gaussian_noise_current_source<RungeKuttaIntegrator>();
}

TEST(noise_source, gaussian_noise_current_source_dormand_prince)
{
	test_gaussian_noise_current_source<DormandPrinceIntegrator>();
}
}
