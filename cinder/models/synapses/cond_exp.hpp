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
 * @file cond_exp.hpp
 *
 * Implementation of a conductance based synapse with exponential decay.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_SYNAPSES_COND_EXP_HPP
#define CINDER_MODELS_SYNAPSES_COND_EXP_HPP

#include <cinder/models/synapse.hpp>

namespace cinder {
/**
 * State vector used for the CondExp synapse. The g_syn() state component holds
 * the conductivity of the synapse.
 */
struct CondExpState : public VectorBase<CondExpState, Real, 1> {
	using VectorBase<CondExpState, Real, 1>::VectorBase;

	TYPED_VECTOR_ELEMENT(g_syn, 0, Conductance);
};

/**
 * Parameter vector used for in the CondExp synapse. A CondExp synapse is
 * described by three parameters. The synaptic weight w_syn() determines the
 * the value by which the conductivity of the synapse is increased when an
 * external spike arrives. The time constant tau_syn() is the time constant
 * of the exponential decay. The e_rev_syn() parameter controls the reversal
 * potential of the ion channel controlled by the synapse.
 */
struct CondExpParameters : public VectorBase<CondExpParameters, Real, 3> {
	using VectorBase<CondExpParameters, Real, 3>::VectorBase;

	/**
	 * Default constructor of the CondExpParameters class. Sets the synaptic
	 * weight to 0.1uS, the synaptic time constant to 5ms and the reversal
	 * potential to 0V.
	 */
	CondExpParameters()
	{
		w_syn(0.1_uS);
		tau_syn(5_ms);
		e_rev_syn(0_mV);
	}

	TYPED_VECTOR_ELEMENT(w_syn, 0, Conductance);
	TYPED_VECTOR_ELEMENT(tau_syn, 1, RealTime);
	TYPED_VECTOR_ELEMENT(e_rev_syn, 2, Voltage);
};

/**
 * The CondExp synapse implements a conductivity based synapse with exponential
 * decay.
 */
struct CondExp : public ConductanceBasedSynapseBase<CondExp, CondExpState,
                                                    CondExpParameters> {
	using Base =
	    ConductanceBasedSynapseBase<CondExp, CondExpState, CondExpParameters>;
	using Base::Base;
	using Base::p;

	friend SynapseBase<CondExp, CondExpState, CondExpParameters>;

private:
	Real m_tau_inv;

	/**
	 * Handles the arival of an input spike -- simply increases the conductivity
	 * of the synapse by increasing the first state component.
	 */
	template <typename State, typename System>
	void process_spike(const Spike &spike, Time, State &s, System &) const
	{
		s[0] += p().w_syn() * spike.w;
	}

public:
	/**
	 * Constructor of the CondExp synapse, allowing to directly set the synapse
	 * parameters.
	 *
	 * @param w_syn is the synaptic weight, determining by how much the
	 * conductance is increased whenever an external spike is received.
	 * @param tau_syn is the synaptic time constant, determining the slope of
	 * the exponential decay.
	 * @param e_rev_syn is the reversal potential of the channels controlled
	 * be the synapse, determining the potential to which the membrane is pulled
	 * whenever the synapse receives input spikes.
	 * @param input_spikes is a list containing the input spike trains.
	 */
	CondExp(Conductance w_syn, RealTime tau_syn, Voltage e_rev_syn,
	        const std::vector<Spike> &input_spikes = std::vector<Spike>())
	    : Base({{w_syn, tau_syn, e_rev_syn}}, input_spikes)
	{
	}

	/**
	 * Calculates the inverse of the tau_syn parameter to speed up the
	 * time-critical calls to df().
	 */
	template <typename State2, typename System>
	void init(Time, const State2 &, const System &)
	{
		m_tau_inv = 1.0_R / p().tau_syn();
	}

	/**
	 * Implements the dynamics of the CondExp synapse which are merely an
	 * exponential decay of the conductivity state component.
	 */
	template <typename State2, typename System>
	State df(const State2 &s, const System &) const
	{
		return State({-m_tau_inv * s[0]});
	}
};
}

#endif /* CINDER_MODELS_SYNAPSES_COND_EXP_HPP */

