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
 * @file cur_exp.hpp
 *
 * Implementatino of a current base synapse with exponential decay.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_SYNAPSES_CUR_EXP_HPP
#define CINDER_MODELS_SYNAPSES_CUR_EXP_HPP

#include <cinder/models/synapse.hpp>

namespace cinder {
/**
 * State vector used by the CurExp synape type.
 */
struct CurExpState : public VectorBase<CurExpState, Real, 1> {
	using VectorBase<CurExpState, Real, 1>::VectorBase;

	TYPED_VECTOR_ELEMENT(i_syn, 0, Current);
};

/**
 * Current based synapse with exponential decay.
 */
struct CurExp : public SynapseBase<CurExp, CurExpState> {
private:
	friend SynapseBase<CurExp, CurExpState>;

	Current m_w;
	Real m_tau_inv;

	template <typename State2, typename System>
	void process_spike(const Spike &spike, Time, State2 &s, System &) const
	{
		// Increase the channel current
		s[0] += m_w * spike.w;
	}

public:
	CurExp(Current w, Time tau,
	       const std::vector<Spike> &input_spikes = std::vector<Spike>())
	    : SynapseBase<CurExp, CurExpState>(input_spikes),
	      m_w(w),
	      m_tau_inv(1.0_R / tau.sec())
	{
	}

	template <typename State2, typename System>
	State df(const State2 &s, const System &) const
	{
		// Exponential decay of the current
		return State({-m_tau_inv * s[0]});
	}
};
}

#endif /* CINDER_MODELS_SYNAPSES_CUR_EXP_HPP */
