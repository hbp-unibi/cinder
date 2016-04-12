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

/**
 * @file solver.hpp
 *
 * Contains the actual class which is responsible for solving autonomous
 * ordinary differential equations which possess discontinuities at certain
 * time points.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_ODE_SOLVER_HPP
#define CINDER_ODE_SOLVER_HPP

#include <tuple>

#include <cinder/common/time.hpp>
#include <cinder/ode/controller.hpp>

namespace cinder {
/**
 * Helper class which orchestrates the integrator, recorder and controller in
 * order to solve the given differential equation.
 */
template <typename ODE, typename Integrator, typename Recorder,
          typename Controller>
class Solver {
public:
	using State = typename ODE::State;

private:
	/**
	 * Compound object which stores all information about the simulation. A
	 * reference to the System object is passed to most ODE member functions for
	 * access to the global state.
	 */
	class System {
	private:
		ODE &m_ode;
		Integrator &m_integrator;
		Recorder &m_recorder;
		Controller &m_controller;
		Time m_t;
		State m_s;

	public:
		System(ODE &ode, Integrator &integrator, Recorder &recorder,
		       Controller &controller, Time t0)
		    : m_ode(ode),
		      m_integrator(integrator),
		      m_recorder(recorder),
		      m_controller(controller),
		      m_t(t0),
		      m_s(ode.s0())
		{
		}

		const ODE &ode() const { return m_ode; }
		ODE &ode() { return m_ode; }

		const Integrator &integrator() const { return m_integrator; }
		Integrator &integrator() { return m_integrator; }

		const Recorder &recorder() const { return m_recorder; }
		Recorder &recorder() { return m_recorder; }

		const Controller &controller() const { return m_controller; }
		Controller &controller() { return m_controller; }

		const State &s() const { return m_s; }
		State &s() { return m_s; }

		Time t() const { return m_t; }
		Time &t() { return m_t; }
		void t(Time t) { m_t = t; }
	};

	System m_sys;

public:
	/**
	 * Constructor of the solver class.
	 */
	Solver(ODE &ode, Integrator &integrator, Recorder &recorder,
	       Controller &controller, Time t0 = 0_s)
	    : m_sys(ode, integrator, recorder, controller, t0)
	{
	}

	/**
	 * Returns the current time the solver is at.
	 */
	Time t() const { return m_sys.t(); }

	/**
	 * Sets the current time the solver is at. The given time must be positive.
	 */
	void t(Time t) { m_sys.t() = t; }

	/**
	 * Returns the state the solver currently is at.
	 */
	const State &s() const { return m_sys.s(); }

	/**
	 * Returns the system instance, allowing access to all underlying types.
	 */
	const System &sys() const { return m_sys; }
	System &sys() { return m_sys; }

	/**
	 * Advances the simulation up to time tEnd with the given time step. Calles
	 * the callback functions in the given Recorder, Integrator and Controller
	 * instances.
	 */
	void solve(Time tEnd = MAX_TIME, Time tDelta = 0.1_ms)
	{
		// Fetch some references at the simulation object for convenience
		Time &t = m_sys.t();
		State &s = m_sys.s();
		ODE &ode = m_sys.ode();
		Recorder &recorder = m_sys.recorder();

		// Iterate over all time slices
		Time dt;
		while (t < tEnd && t >= Time(0)) {
			// Fetch the time of the next discontinuity and integrate at most up
			// to this time
			const Time tNextDiscontinuity =
			    std::max(t, std::min(tEnd, const_cast<const ODE &>(ode)
			                                   .next_discontinuity(t)));
			if (tNextDiscontinuity != t) {
				const Time tDeltaMax = tNextDiscontinuity - t;
				std::tie(s, dt) = m_sys.integrator().integrate(
				    std::min(tDelta, tDeltaMax), tDeltaMax,
				    const_cast<const State &>(s),
				    const_cast<const System &>(m_sys));
				t += dt;

				// Call the system update function
				ode.update(const_cast<const Time &>(t), s, m_sys);

				// Record the updated state
				recorder.record(const_cast<const Time &>(t),
				                const_cast<const State &>(s),
				                const_cast<const System &>(m_sys));
			}

			// Once the discontinuity is reached, call the ODE update function,
			// record the state again in order to properly represent the sharp
			// jump in the state vector
			if (t >= tNextDiscontinuity) {
				ode.handle_discontinuity(const_cast<const Time &>(t), s, m_sys);
				t += Time(1);  // Advance by the smallest possible time-step
				recorder.record(const_cast<const Time &>(t),
				                const_cast<const State &>(s),
				                const_cast<const System &>(m_sys));
			}

			// Ask the controller whether it is time to abort
			const ControllerResult cres = m_sys.controller().control(
			    const_cast<const Time &>(t), const_cast<const State &>(s),
			    const_cast<const System &>(m_sys));
			if (cres == ControllerResult::ABORT ||
			    (cres == ControllerResult::MAY_CONTINUE &&
			     tNextDiscontinuity == tEnd)) {
				return;
			}
		}
	}
};

/**
 * Convenience function for generating a solver for the given ODE, Integrator,
 * Recorder and Controller.
 */
template <typename ODE, typename Integrator, typename Recorder,
          typename Controller>
static Solver<ODE, Integrator, Recorder, Controller> make_solver(
    ODE &ode, Integrator &integrator, Recorder &recorder,
    Controller &controller, Time t0 = 0_s)
{
	return Solver<ODE, Integrator, Recorder, Controller>(
	    ode, integrator, recorder, controller, t0);
}
}

#endif /* CINDER_ODE_SOLVER_HPP */
