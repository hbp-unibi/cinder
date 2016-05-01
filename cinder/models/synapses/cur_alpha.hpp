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
struct AlphaCurrentState : public VectorBase<AlphaCurrentState, Real, 2> {
	using VectorBase<AlphaCurrentState, Real, 2>::VectorBase;

	static constexpr AlphaCurrentState scale()
	{
		return AlphaCurrentState({1e9, 1e9});
	}
};

/**
 * Alpha-function shaped current-based synapse.
 */
struct CurAlpha : public SynapseBase<CurAlpha, AlphaCurrentState> {
private:
	friend SynapseBase<CurAlpha, AlphaCurrentState>;

	Real m_w;
	Real m_tau_inv;

	template <typename State2, typename System>
	void process_spike(const Spike &spike, Time, State2 &s, System &) const
	{
		s[1] += m_w * spike.w;
	}

public:
	CurAlpha(Current w, Time tau,
	         const std::vector<Spike> &input_spikes = std::vector<Spike>())
	    : SynapseBase<CurAlpha, AlphaCurrentState>(input_spikes),
	      m_w(w.v() * 2.718281828_R),
	      m_tau_inv(1.0_R / tau.sec())
	{
	}

	template <typename State2, typename System>
	State df(const State2 &s, const System &) const
	{
		return State({(s[1] - s[0]) * m_tau_inv, -s[1] * m_tau_inv});
	}
};
}

#endif /* CINDER_MODELS_SYNAPSES_CUR_ALPHA_HPP */
