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

#include <cinder/common/time.hpp>
#include <cinder/common/types.hpp>

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
 * Controller type which runs the simulation until the neuron settles to its
 * resting potential without any currents flowing.
 */
struct NeuronController {
	template <typename State, typename System>
	static ControllerResult control(Time, const State &s, const System &sys)
	{
		static constexpr Voltage MAX_DELTA_V = 1_uV;
		static constexpr Current MAX_DELTA_I = 1_pA;

		// Abort if there are no more input spikes, the neuron membrane voltage
		// is near the resting potential and the current is near zero.
		if (std::abs(sys.ode().voltage(s, sys) -
		             sys.ode().membrane().p().v_rest()) < MAX_DELTA_V.v() &&
		    std::abs(sys.ode().current(s, sys) < MAX_DELTA_I.v())) {
			return ControllerResult::MAY_CONTINUE;
		}
		return ControllerResult::CONTINUE;  // Go on forever
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
 * Class used to cascade a number of Recorders. Use the makeMultiRecorder()
 * method to conveniently construct a Recorder consisting of multiple recorders.
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
}

#endif /* CINDER_ODE_CONTROLLER_HPP */
