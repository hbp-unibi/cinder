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
 * @file euler.hpp
 *
 * Contains an implementation of the basic Euler integrator.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_INTEGRATOR_EULER_HPP
#define CINDER_INTEGRATOR_EULER_HPP

#include <utility>

#include <cinder/common/time.hpp>
#include <cinder/common/types.hpp>

namespace cinder {
/**
 * The Euler class represents euler's method for integrating ODEs.
 */
class EulerIntegrator {
public:
	/**
	 * Implements the Euler integration method.
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
		return std::make_pair(s + h * sys.ode().df(s, sys), tDelta);
	}
};
}

#endif /* CINDER_INTEGRATOR_EULER_HPP */

