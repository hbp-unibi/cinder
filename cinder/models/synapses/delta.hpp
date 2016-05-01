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
 * @file delta.hpp
 *
 * Synapse which implements a simple dirac-delta current pulse (voltage step).
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_SYNAPSES_DELTA_HPP
#define CINDER_MODELS_SYNAPSES_DELTA_HPP

#include <cinder/models/synapse.hpp>

namespace cinder {
/**
 * Synapse type which increases the membrane potential by a fixed value.
 */
struct Delta : public SynapseBase<Delta, NullState> {
private:
	friend SynapseBase<Delta, NullState>;

	Voltage m_pulse_v;

	template <typename State, typename System>
	void process_spike(const Spike &spike, Time, State &, System &sys) const
	{
		// The first system state is supposed to be the membrane voltage
		sys.s()[0] += m_pulse_v * spike.w;
	}

public:
	Delta(Voltage pulse_v,
	        const std::vector<Spike> &input_spikes = std::vector<Spike>())
	    : SynapseBase<Delta, NullState>(input_spikes),
	      m_pulse_v(pulse_v)
	{
	}

	template <typename State2, typename System>
	State df(const State2 &, const System &) const
	{
		return NullState();
	}

	template <typename State, typename System>
	Current current(const State &, const System &) const
	{
		return Current();
	}
};
}

#endif /* CINDER_MODELS_SYNAPSES_DELTA_HPP */
