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
 * @file midpoint.hpp
 *
 * Contains an implementation of the midpoint method.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_INTEGRATOR_MIDPOINT_HPP
#define CINDER_INTEGRATOR_MIDPOINT_HPP

#include <utility>

#include <cinder/common/time.hpp>
#include <cinder/common/types.hpp>

namespace cinder {
/**
 * The Midpoint class implements the second-order Runge-Kutta method.
 */
class MidpointIntegrator {
public:
	/**
	 * Implements the second-order Runge-Kutta method (Midpoint method).
	 *
	 * @param tDelta is the timestep width.
	 * @param s is the current state vector at the previous timestep.
	 * @param df is the function which calculates the derivative for a given
	 * state.
	 * @return the new state for the next timestep and the actually used
	 * timestep.
	 */
	template <typename State, typename System>
	static std::pair<State, Time> integrate(Time tDelta, Time, const State &s,
	                                        const System &sys)
	{
		const Real h = tDelta.sec();
		const State k1 = h * sys.ode().df(s, sys);
		const State k2 = h * sys.ode().df(s + 0.5f * k1, sys);
		return std::pair<State, Time>(s + k2, tDelta);
	}
};
}

#endif /* CINDER_INTEGRATOR_MIDPOINT_HPP */

