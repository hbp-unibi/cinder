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
 * @file cond_alpha.hpp
 *
 * Implementation of a conductance based synapse with alpha-function shaped
 * conductance.
 *
 * The alpha function has the following form:
 *
 *      a(t) = alpha / tau * t * exp(1 - t / tau)
 *
 * This function can be rewritten as the following system of differential
 * equations:
 *
 *      d/dt b(t) = - b(t) / tau
 *      d/dt a(t) = (a(t) - b(t)) / tau
 *
 * Where the first differential equation is equivalent to an exponential-based
 * synapse and the second equation is a first-order low-pass filter.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_SYNAPSES_COND_ALPHA_HPP
#define CINDER_MODELS_SYNAPSES_COND_ALPHA_HPP

#include <cinder/common/fast_math.hpp>
#include <cinder/models/synapse.hpp>

namespace cinder {
/**
 * State vector used in the CondAlpha synapse. The first state component holds
 * the low-pass filtered conductance, the second state component the conductance
 * without low-pass filter.
 */
struct CondAlphaState : public VectorBase<CondAlphaState, Real, 2> {
	using VectorBase<CondAlphaState, Real, 2>::VectorBase;

	TYPED_VECTOR_ELEMENT(g_syn, 0, Conductance);
	TYPED_VECTOR_ELEMENT(g_syn_exp, 1, Conductance);
};

/**
 * Parameter vector used for the CondAlpha synapse. A CondAlpha synapse posseses
 * three parameters. The synaptic weight w_syn() determines the peak
 * conductivity of the synapse when an external spike arrives. The time constant
 * tau_syn() controls how long the synapse takes to reach its maximum
 * conductivity and the slope of the decay. The e_rev_syn() parameter controls
 * the reversal potential of the ion channel controlled by the synapse.
 */
struct CondAlphaParameters : public VectorBase<CondAlphaParameters, Real, 3> {
	using VectorBase<CondAlphaParameters, Real, 3>::VectorBase;

	/**
	 * Default constructor of the CondAlphaParameters class. Sets the synaptic
	 * weight to 0.1uS, the synaptic time constant to 5ms and the reversal
	 * potential to 0V.
	 */
	CondAlphaParameters()
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
 * Alpha-function shaped conductance-based synapse.
 */
struct CondAlpha : public ConductanceBasedSynapseBase<CondAlpha, CondAlphaState,
                                                      CondAlphaParameters> {
	using Base = ConductanceBasedSynapseBase<CondAlpha, CondAlphaState,
	                                         CondAlphaParameters>;
	using Base::Base;
	using Base::p;

	friend SynapseBase<CondAlpha, CondAlphaState,
	                                   CondAlphaParameters>;

private:
	Real m_tau_inv;

	template <typename State2, typename System>
	void process_spike(const Time &, Time, State2 &s, System &) const
	{
		s[1] += p().w_syn() * 2.718281828_R;
	}

public:
	/**
	 * Constructor allowing to directly set all three synapse parameters.
	 *
	 * @param w_syn is the peak conductivity that is reached when the synapse
	 * receives a single input spike.
	 * @param tau_syn is the synaptic time constant, determining the slope of
	 * the conductivity rise and fall and the relative position of the maximum
	 * conductivity.
	 * @param e_rev_syn is the reversal potential of the channels controlled
	 * be the synapse, determining the potential to which the membrane is pulled
	 * whenever the synapse receives input spikes.
	 * @param input_spikes is a list containing the input spike trains.
	 */
	CondAlpha(Conductance w_syn, RealTime tau_syn, Voltage e_rev_syn,
	          const std::vector<Time> &input_spikes = std::vector<Time>())
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
	 * The df() method describes the dynamics of the CondAlpha synapse,
	 * implementing a low-pass filtered version of the second state component
	 * for the first state component and an exponential decay of the second
	 * state component:
	 *
	 *      d/dt b(t) = - b(t) / tau
	 *      d/dt a(t) = (a(t) - b(t)) / tau
	 */
	template <typename State2, typename System>
	State df(const State2 &s, const System &) const
	{
		return State({(s[1] - s[0]) * m_tau_inv, -s[1] * m_tau_inv});
	}
};
}

#endif /* CINDER_MODELS_SYNAPSES_COND_ALPHA_HPP */

