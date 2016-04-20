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
 * @file ode.hpp
 *
 * Contains a basic implementation of the ODE concept.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_ODE_ODE_HPP
#define CINDER_ODE_ODE_HPP

#include <cinder/common/time.hpp>

namespace cinder {
/**
 * The ODEBase represents a basic ordinary differential equation with no
 * discontinuities.
 *
 * @tparam State is the state vector underlying the ODE.
 */
template <typename State_>
struct ODEBase {
	using State = State_;

	/**
	 * Allows the implementations to perform some initialization.
	 */
	template <typename State, typename System>
	static void init(Time, const State &, const System &)
	{
	}

	/**
	 * Returns the initial state of the current source.
	 */
	static State s0() { return State(); }

	/**
	 * Returns the position of the next discontinuity in the differential
	 * equation.
	 */
	static Time next_discontinuity(Time) { return MAX_TIME; }

	/**
	 * Called whenever the discontinuity is reached.
	 */
	template <typename State, typename System>
	void handle_discontinuity(Time, const State &, const System &)
	{
	}

	/**
	 * Gives the model the opportunity to perform some additional calculations.
	 */
	template <typename State, typename System>
	static void update(Time, const State &, const System &)
	{
	}

	/**
	 * Returns the differential.
	 */
	template <typename State, typename System>
	static State_ df(const State &, const System &)
	{
		return State_();
	}
};
}

#endif /* CINDER_ODE_ODE_HPP */
