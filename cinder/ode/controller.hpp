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
 * Controller type which does not prematurely abort the neuron simulation.
 */
struct NullController {
	template <typename State, typename System>
	static ControllerResult control(Time, const State &, const System &)
	{
		return ControllerResult::CONTINUE;  // Go on forever
	}
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
}

#endif /* CINDER_ODE_CONTROLLER_HPP */
