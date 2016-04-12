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
 * @file current_source.hpp
 *
 * Implementation of some basic current sources. Furthermore implements the
 * MultiCurrentSource object which allows to combine an arbitrary number of
 * current sources into a single current source.
 *
 * A current source in itself is an ODE object which can be integrated using
 * the Solver class.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_SYNAPSE_HPP
#define CINDER_MODELS_SYNAPSE_HPP

#include <queue>

#include <cinder/common/spike.hpp>
#include <cinder/models/current_source.hpp>

namespace cinder {
/**
 * Base class from which synapse implementations should derive. Implements the
 * handling
 */
template <typename Impl, typename StateImpl>
struct SynapseBase : public CurrentSourceBase<StateImpl> {
private:
	std::priority_queue<Spike, std::vector<Spike>, std::greater<Spike>>
	    m_input_spikes;

	Impl &impl() { return static_cast<Impl &>(*this); }

public:
	SynapseBase(const std::vector<Spike> &input_spikes = std::vector<Spike>())
	    : m_input_spikes(input_spikes.begin(), input_spikes.end())
	{
	}

	auto &input_spikes() { return m_input_spikes; }
	const auto &input_spikes() const { return m_input_spikes; }

	/**
	 * Returns the position of the next discontinuity in the differential
	 * equation.
	 */
	Time next_discontinuity(Time) const
	{
		return m_input_spikes.empty() ? MAX_TIME : m_input_spikes.top().t;
	}

	/**
	 * Called whenever the discontinuity is reached. Relays the call to the
	 * child process_spike() method.
	 */
	template <typename State, typename System>
	void handle_discontinuity(Time t, State &s, System &sys)
	{
		while (!m_input_spikes.empty() && m_input_spikes.top().t <= t) {
			impl().process_spike(m_input_spikes.top(), t, s, sys);
			m_input_spikes.pop();
		}
	}
};
}

#endif /* CINDER_MODELS_SYNAPSE_HPP */
