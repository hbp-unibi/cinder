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
 * @file controller.hpp
 *
 * Contains some standard controllers -- the class which is responsible for
 * prematurely aborting the simulation if nothing happens.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_ODE_CONTROLLER_HPP
#define CINDER_ODE_CONTROLLER_HPP

#include <tuple>
#include <vector>

#include <cinder/common/time.hpp>
#include <cinder/common/types.hpp>
#include <cinder/common/vector.hpp>

namespace cinder {
/**
 * Enum describing the descision made by the controller about whether to
 * continue integrating or not.
 */
enum class ControllerResult {
	/**
     * The controller should not abort, unless the end time tEnd is reached.
     */
	CONTINUE,

	/**
     * The controller may abort if there are no more discontinuities in the
     * ODE.
     */
	MAY_CONTINUE,

	/**
     * The controller must abort.
     */
	ABORT
};

/**
 * Controller class which returns a constant result.
 *
 * @tparam Result is the ControllerResult that should be constantly returned
 * by the controller.
 */
template <ControllerResult Result>
struct ConstantController {
	template <typename State, typename System>
	static ControllerResult control(Time, const State &, const System &)
	{
		return Result;
	}
};

/**
 * Controller type which does not prematurely abort the neuron simulation.
 */
struct NullController : public ConstantController<ControllerResult::CONTINUE> {
};

/**
 * Controller type which runs the simulation until the system of differntial
 * equations has converged to a resting state.
 */
class AutoController {
private:
	Time m_last_time;
	std::vector<Real> m_last_state;

public:
	/**
	 * Allows premature abortion of the differential equation integration once
	 * the overall system activity sinks below a certain threshold relative to
	 * the current value range of the system.
	 */
	template <typename State, typename System>
	ControllerResult control(Time t, const State &s, const System &)
	{
		static constexpr Real MAX_ACTIVITY = 1e-3;
		static constexpr Real MAX_ACTIVITY_REL = 1e-3;

		// Calculate the activity of the neuron by numerically calculating the
		// differential
		bool tripped = false;
		if (t - m_last_time > 1_ms) {
			if (!m_last_state.empty()) {
				const Time dt = t - m_last_time;
				const State ds =
				    (s - VectorView<Real, State::Size>(m_last_state.data())) *
				    State::scales() / dt.sec();
				const Real activity = ds.L2Norm();

				// Calculate the current threshold, taking the absolute scale of
				// the state vector into account
				const Real threshold =
				    MAX_ACTIVITY +
				    (s * State::scales()).L2Norm() * MAX_ACTIVITY_REL;

				tripped = activity < threshold;
			}

			// Copy the current state and timestamp to the last_state and
			// last_time variables
			m_last_state.assign(s.begin(), s.end());
			m_last_time = t;
		}

		return tripped ? ControllerResult::MAY_CONTINUE
		               : ControllerResult::CONTINUE;
	}
};

/**
 * Controller class which can be used to describe an external abort condition
 * using a lambda expression.
 *
 * @tparam F is the type of the lambda expression.
 * @tparam DefaultResult is the result that should be returned when the lambda
 * returns true.
 */
template <typename F,
          ControllerResult DefaultResult = ControllerResult::MAY_CONTINUE>
struct ConditionedController {
private:
	F m_f;

public:
	ConditionedController(F f) : m_f(std::move(f)) {}
	template <typename State, typename System>
	ControllerResult control(Time, const State &, const System &)
	{
		if (m_f()) {
			return DefaultResult;
		}
		return ControllerResult::ABORT;
	}
};

/**
 * Creates a new ConditionedController instance.
 */
template <typename F,
          ControllerResult DefaultResult = ControllerResult::MAY_CONTINUE>
static ConditionedController<F, DefaultResult> make_conditioned_controller(
    const F &f)
{
	return ConditionedController<F, DefaultResult>(f);
}

/**
 * Lowest recursion level of the MultiController class. Returns "MAY_CONTINUE".
 */
template <typename... Controllers>
class MultiController
    : public ConstantController<ControllerResult::MAY_CONTINUE> {
};

/**
 * Class used to cascade a number of Controllers. Use the
 * make_multi_controller()
 * method to conveniently construct a Controller consisting of multiple
 * controllers.
 */
template <typename Controller, typename... Controllers>
class MultiController<Controller, Controllers...>
    : MultiController<Controllers...> {
private:
	Controller &controller;

public:
	MultiController(Controller &controller, Controllers &... cs)
	    : MultiController<Controllers...>(cs...), controller(controller)
	{
	}

	template <typename State, typename System>
	ControllerResult control(Time t, const State &s, const System &sys)
	{
		// Evaluate this controller -- abort if it forces abortion
		const ControllerResult res1 = controller.control(t, s, sys);
		if (res1 == ControllerResult::ABORT) {
			return ControllerResult::ABORT;
		}

		// Evaluate the other controllers -- relay abortion
		const ControllerResult res2 =
		    MultiController<Controllers...>::control(t, s, sys);
		if (res2 == ControllerResult::ABORT) {
			return ControllerResult::ABORT;
		}

		// If one controller wants to continue, continue
		if (res1 == ControllerResult::CONTINUE ||
		    res2 == ControllerResult::CONTINUE) {
			return ControllerResult::CONTINUE;
		}

		// Otherwise output "MAY_CONTINUE"
		return ControllerResult::MAY_CONTINUE;
	}
};

/**
 * Helper function which can be used to conveniently create a MultiRecorder
 * instance. Simply pass references to all recorders that should be used to the
 * method and store the result in an auto variable.
 */
template <typename... Controllers>
MultiController<Controllers...> make_multi_controller(Controllers &... cs)
{
	return MultiController<Controllers...>(cs...);
}

/**
 * The ConditionedAutoController class is a combination of the
 * "AutoController" and the "ConditionedController" classes. It allows to
 * simulate a neuron until it has settled or an externally defined condition is
 * reached. Use the make_conditioned_auto_controller() method to conveniently
 * create an instance of the ConditionedAutoController.
 *
 * @tparam F is the
 */
template <typename F>
struct ConditionedAutoController
    : MultiController<ConditionedController<F>, AutoController> {
	ConditionedController<F> conditioned_controller;
	AutoController auto_controller;

	ConditionedAutoController(F f)
	    : MultiController<ConditionedController<F>, AutoController>(
	          conditioned_controller, auto_controller),
	      conditioned_controller(f)
	{
	}
};

/**
 * Returns an instance of the ConditionedAutoController class, which is a
 * combination of a AutoController and a ConditionedController.
 *
 * @param f is the function describing the abort condition.
 */
template <typename F>
ConditionedAutoController<F> make_conditioned_auto_controller(F f)
{
	return ConditionedAutoController<F>(f);
}
}

#endif /* CINDER_ODE_CONTROLLER_HPP */
