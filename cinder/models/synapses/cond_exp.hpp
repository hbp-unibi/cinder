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
 * State vector used for the conductance based neuron.
 */
struct SingleConductanceState
    : public VectorBase<SingleConductanceState, Real, 1> {
	using VectorBase<SingleConductanceState, Real, 1>::VectorBase;

	static constexpr SingleConductanceState scale() { return SingleConductanceState({1e6}); }
};

/**
 * Current based synapse with exponential decay.
 */
struct CondExp : public SynapseBase<CondExp, SingleConductanceState> {
private:
	friend SynapseBase<CondExp, SingleConductanceState>;

	Conductance m_w;
	Real m_tau_inv;
	Voltage m_e_rev;

	template <typename State, typename System>
	void process_spike(const Spike &spike, Time, State &s, System &) const
	{
		// Increase the channel conductance
		s[0] += m_w * spike.w;
	}

public:
	CondExp(Conductance w, Time tau, Voltage e_rev,
	        const std::vector<Spike> &input_spikes = std::vector<Spike>())
	    : SynapseBase<CondExp, SingleConductanceState>(input_spikes),
	      m_w(w),
	      m_tau_inv(1.0_R / tau.sec()),
	      m_e_rev(e_rev)
	{
	}

	template <typename State, typename System>
	SingleConductanceState df(const State &s, const System &) const
	{
		// Exponential decay of the conductance
		return SingleConductanceState({-m_tau_inv * s[0]});
	}

	template <typename State, typename System>
	Current current(const State &s, const System &sys) const
	{
		// Calculate the current depending on the channel conductance and the
		// current voltage
		return Current(s[0] * (m_e_rev - sys.ode().voltage(sys.s(), sys)).v());
	}
};
}

#endif /* CINDER_MODELS_SYNAPSES_COND_EXP_HPP */
