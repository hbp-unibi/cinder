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
 * @file dirac.hpp
 *
 * Synapse which implements a simple current pulse (voltage step).
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_SYNAPSES_DIRAC_HPP
#define CINDER_MODELS_SYNAPSES_DIRAC_HPP

#include <cinder/models/synapse.hpp>

namespace cinder {
/**
 * Synapse type which increases the membrane potential by a fixed value.
 */
struct Dirac : public SynapseBase<Dirac, NullState> {
private:
	friend SynapseBase<Dirac, NullState>;

	Voltage m_pulse_v;

	template <typename State, typename System>
	void process_spike(const Spike &spike, Time, State &, System &sys) const
	{
		// The first system state is supposed to be the membrane voltage
		sys.s()[0] += m_pulse_v * spike.w;
	}

public:
	Dirac(Voltage pulse_v,
	        const std::vector<Spike> &input_spikes = std::vector<Spike>())
	    : SynapseBase<Dirac, NullState>(input_spikes),
	      m_pulse_v(pulse_v)
	{
	}

	template <typename State, typename System>
	NullState df(const State &, const System &) const
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

#endif /* CINDER_MODELS_SYNAPSES_DIRAC_HPP */
