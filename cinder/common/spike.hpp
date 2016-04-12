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
 * @file spike.hpp
 *
 * Contains the data type describing a spike and auxiliary functions for
 * generating and scaling input spikes.
 *
 * @author Andreas Stöckel
 */

#ifndef CINDER_COMMON_SPIKE_HPP
#define CINDER_COMMON_SPIKE_HPP

#include <cinder/common/types.hpp>
#include <cinder/common/time.hpp>

namespace cinder {
/**
 * Structure representing a single input spike.
 */
struct Spike {
	/**
	 * Enum describing the two possible spike types: Inhibitory and excitatory
	 * spikes. Generally spike which arrive at a synapse with a positive
	 * synaptic weight are excitatory, spikes which arrive at a synapse with
	 * negative synaptic weight are inhibitory.
	 */
	enum class Type {
		INHIBITORY, EXCITATORY
	};

	/**
	 * Time at which the spike is received by the neuron.
	 */
	Time t;

	/**
	 * Weight of the spike -- can be understood as the weight of the synaptic
	 * connection the spike arrives at.
	 */
	Real w;

	/**
	 * Constructor of the Spike class, allows to initialize all members.
	 *
	 * @param t is the time at which the spike is received by the neuron.
	 * @param w is the weight of the spike.
	 */
	constexpr Spike(Time t = Time(0), Real w = 1.0) : t(t), w(w) {}

	/**
	 * Returns whether this is an inhibitory or excitatory spike.
	 */
	Type type() const {
		return (w < 0.0) ? Type::INHIBITORY : Type::EXCITATORY;
	}

	/**
	 * Operator used to sort spikes by time. This way a vector of spikes can be
	 * easily sorted using the STL methods.
	 */
	friend bool operator<(const Spike &s1, const Spike &s2)
	{
		return s1.t < s2.t;
	}

	/**
	 * Operator used to sort spikes by time. This way a vector of spikes can be
	 * easily sorted using the STL methods.
	 */
	friend bool operator>(const Spike &s1, const Spike &s2)
	{
		return s1.t > s2.t;
	}
};
}

#endif /* CINDER_COMMON_SPIKE_HPP */
