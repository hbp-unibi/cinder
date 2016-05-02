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
 * Parameters of a Delta synapse. Determines by how much the membrane potential
 * is increased/decreased whenever an input spike is received.
 */
struct DeltaParameters : public VectorBase<DeltaParameters, Real, 1> {
	using VectorBase<DeltaParameters, Real, 1>::VectorBase;

	TYPED_VECTOR_ELEMENT(w_syn, 0, Voltage);
};

/**
 * Synapse type which increases the membrane potential by a fixed value, which
 * corresponds to the Dirac-Delta shaped current pulse (hence the name of the
 * synapse).
 */
struct Delta : public SynapseBase<Delta, NullState, DeltaParameters> {
	using Base = SynapseBase<Delta, NullState, DeltaParameters>;
	using Base::Base;
	using Base::p;

	friend SynapseBase<Delta, NullState, DeltaParameters>;

private:
	template <typename State, typename System>
	void process_spike(const Spike &spike, Time, State &, System &sys) const
	{
		// The first system state is supposed to be the membrane voltage
		sys.s()[0] += p().w_syn() * spike.w;
	}

public:
	/**
	 * Constructor of the Delta synapse, allowing to directly set the synapse
	 * parameters.
	 *
	 * @param w_syn is the synaptic weight, determining the voltage by which
	 * the membrane potential is increased/decreased whenever a spike arrives
	 * at the synapse.
	 * @param input_spikes is a list containing the input spike trains.
	 */
	Delta(Voltage w_syn,
	        const std::vector<Spike> &input_spikes = std::vector<Spike>())
	    : Base({{w_syn}}, input_spikes)
	{
	}

	/**
	 * As the Delta synapse does not induce a current (appart from a dirac-delta
	 * shaped current, which is not modelled here) and does not possess a state
	 * vector, the "current" method is overriden to simply return zero.
	 */
	template <typename State2, typename System>
	static Current current(const State2 &, const System &)
	{
		return Current();
	}
};
}

#endif /* CINDER_MODELS_SYNAPSES_DELTA_HPP */

