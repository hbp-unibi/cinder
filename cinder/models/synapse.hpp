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
 * @file synapse.hpp
 *
 * Contains the base class of all synapses, SynapseBase. A synapse is a current
 * source modulated by a set of input spikes.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_SYNAPSE_HPP
#define CINDER_MODELS_SYNAPSE_HPP

#include <queue>

#include <cinder/models/current_source.hpp>

namespace cinder {
/**
 * Base class from which synapse implementations should derive. Implements the
 * handling of input spikes.
 *
 * @tparam Impl_ is the top-level class which actually implements the synaptic
 * dynamics and possesses a process_spike() method, which is called whenever an
 * input spike is received.
 * @tparam State_ is a vector type describing the current synapse state in its
 * corresponding differential equation.
 * @tparam Parameters_ is a vector type holding the user settable synaptic
 * parameters.
 */
template <typename Impl_, typename State_, typename Parameters_>
struct SynapseBase : public CurrentSourceBase<State_, Parameters_> {
private:
	/**
	 * Priority queue holding the input spikes, with the earliest spikes being
	 * at the beginnign of the queue.
	 */
	std::priority_queue<Time, std::vector<Time>, std::greater<Time>>
	    m_input_spikes;

public:
	using Base = CurrentSourceBase<State_, Parameters_>;

	using Parameters = Parameters_;
	using State = State_;

	/**
	 * Constructor of the SynapseBase class, copies the given parameters and
	 * input spikes into internal buffers.
	 *
	 * @param params is an instance of the Parameters vector from which the
	 * user-defined synaptic parameters are read.
	 * @param input_spikes is a list containing the input spikes. The spikes are
	 * automatically sorted into a priority queue and may be unsorted.
	 */
	SynapseBase(const Parameters &params = Parameters(),
	            const std::vector<Time> &input_spikes = std::vector<Time>())
	    : Base(params), m_input_spikes(input_spikes.begin(), input_spikes.end())
	{
	}

	/**
	 * Returns a reference at the priority queue holding the input spikes.
	 */
	auto &input_spikes() { return m_input_spikes; }

	/**
	 * Returns a const-reference at the priority queue holding the input spikes.
	 */
	const auto &input_spikes() const { return m_input_spikes; }

	/**
	 * Returns the position of the next discontinuity in the differential
	 * equation.
	 */
	Time next_discontinuity(Time) const
	{
		return m_input_spikes.empty() ? MAX_TIME : m_input_spikes.top();
	}

	/**
	 * Called whenever the discontinuity is reached. Relays the call to the
	 * child process_spike() method.
	 */
	template <typename State2, typename System>
	void handle_discontinuity(Time t, State2 &s, System &sys)
	{
		while (!m_input_spikes.empty() && m_input_spikes.top() <= t) {
			static_cast<Impl_ &>(*this)
			    .process_spike(m_input_spikes.top(), t, s, sys);
			m_input_spikes.pop();
		}
	}
};

/**
 * Base class for conductance based synapses. Contains a current function which
 * converts the conductance stored in the first state component into a current
 * based on the current membrane potential.
 */
template <typename Impl_, typename State_, typename Parameters_>
struct ConductanceBasedSynapseBase
    : public SynapseBase<Impl_, State_, Parameters_> {
	using Base = SynapseBase<Impl_, State_, Parameters_>;
	using Base::Base;
	using Base::p;

	using Parameters = Parameters_;
	using State = State_;

	/**
	 * Converts the conductance stored in the first state component into a
	 * current depending on the current membrane potential and the e_rev_syn()
	 * parameter.
	 */
	template <typename State2, typename System>
	Current current(const State2 &s, const System &sys) const
	{
		// Calculate the current depending on the channel conductance and the
		// current voltage
		return Current(s[0] *
		               (p().e_rev_syn() - sys.ode().voltage(sys.s(), sys)).v());
	}
};

/**
 * Base class for current based synapses. Assumes that the current is stored
 * in the first state component (as is the default for all classes derived from
 * CurrentSourceBase). This class does nothing special, but is provided for
 * symmetry reasons along with ConducdanceBasedSynapseBase.
 */
template <typename Impl_, typename State_, typename Parameters_>
struct CurrentBasedSynapseBase
    : public SynapseBase<Impl_, State_, Parameters_> {
	using Base = SynapseBase<Impl_, State_, Parameters_>;
	using Base::Base;
};
}

#endif /* CINDER_MODELS_SYNAPSE_HPP */

