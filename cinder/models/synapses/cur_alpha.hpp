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
 * @file cur_alpha.hpp
 *
 * Implementation of a synapse which injects an alpha-function shaped current
 * into the neuron.
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

#ifndef CINDER_MODELS_SYNAPSES_CUR_ALPHA_HPP
#define CINDER_MODELS_SYNAPSES_CUR_ALPHA_HPP

#include <cinder/common/fast_math.hpp>
#include <cinder/models/synapse.hpp>

namespace cinder {
/**
 * State vector used in the CurAlpha synapse. The first state component holds
 * the low-pass filtered current, the second state component the current without
 * low-pass filter.
 */
struct CurAlphaState : public VectorBase<CurAlphaState, Real, 2> {
	using VectorBase<CurAlphaState, Real, 2>::VectorBase;

	TYPED_VECTOR_ELEMENT(i_syn, 0, Current);
	TYPED_VECTOR_ELEMENT(i_syn_exp, 1, Current);
};

/**
 * Parameter vector used for the CurAlpha synapse. A CurAlpha synapse posseses
 * two parameters. The synaptic weight w_syn() determines the peak
 * current of the synapse as a result to an external spike The time constant
 * tau_syn() controls how long the synapse takes to reach its maximum
 * conductivity, as well as the slope of the decay.
 */
struct CurAlphaParameters : public VectorBase<CurAlphaParameters, Real, 2> {
	using VectorBase<CurAlphaParameters, Real, 2>::VectorBase;

	/**
	 * Default constructor of the CurAlphaParameters class. Sets the synaptic
	 * weight to 0.1nA and the synaptic time constant to 5ms.
	 */
	CurAlphaParameters()
	{
		w_syn(0.1_nA);
		tau_syn(5_ms);
	}

	TYPED_VECTOR_ELEMENT(w_syn, 0, Current);
	TYPED_VECTOR_ELEMENT(tau_syn, 1, RealTime);
};

/**
 * Alpha-function shaped current-based synapse.
 */
struct CurAlpha : public CurrentBasedSynapseBase<CurAlpha, CurAlphaState,
                                                 CurAlphaParameters> {
	using Base = CurrentBasedSynapseBase<CurAlpha, CurAlphaState,
	                                         CurAlphaParameters>;
	using Base::Base;
	using Base::p;

	friend SynapseBase<CurAlpha, CurAlphaState, CurAlphaParameters>;

private:
	Real m_tau_inv;

	template <typename State2, typename System>
	void process_spike(const Spike &spike, Time, State2 &s, System &) const
	{
		s[1] += p().w_syn() * spike.w * 2.718281828_R;
	}

public:
	/**
	 * Constructor of the CurAlpha synapse, allowing to directly set the synapse
	 * parameters.
	 *
	 * @param w_syn is the synaptic weight, determining the current peak when
	 * a single external spike arrives.
	 * @param tau_syn is the time constant, determining the slope of the current
	 * rise and fall.
	 * @param input_spikes is a list containing the input spike trains.
	 */
	CurAlpha(Current w_syn, RealTime tau_syn,
	         const std::vector<Spike> &input_spikes = std::vector<Spike>())
	    : Base({{w_syn, tau_syn}}, input_spikes)
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
	 * The df() method describes the dynamics of the CurAlpha synapse,
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

#endif /* CINDER_MODELS_SYNAPSES_CUR_ALPHA_HPP */

