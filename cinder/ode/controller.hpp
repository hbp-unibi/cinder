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
}

#endif /* CINDER_ODE_CONTROLLER_HPP */
