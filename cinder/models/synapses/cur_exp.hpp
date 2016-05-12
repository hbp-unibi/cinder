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
 * State vector used by the CurExp synape type, stores the current induced by
 * the synapse in the first state variable.
 */
struct CurExpState : public VectorBase<CurExpState, Real, 1> {
	using VectorBase<CurExpState, Real, 1>::VectorBase;

	TYPED_VECTOR_ELEMENT(i_syn, 0, Current);
};

/**
 * Parameter vector used for the CurExp synapse. A CurExp synapse posseses
 * two parameters. The synaptic weight w_syn() determines by how much the
 * current induced by the synapse is increased whenever an external spike
 * arrives. The parameter tau_syn() is time constant of the exponential decay.
 */
struct CurExpParameters : public VectorBase<CurExpParameters, Real, 2> {
	using VectorBase<CurExpParameters, Real, 2>::VectorBase;

	CurExpParameters()
	{
		w_syn(10_nA);
		tau_syn(5_ms);
	}

	TYPED_VECTOR_ELEMENT(w_syn, 0, Current);
	TYPED_VECTOR_ELEMENT(tau_syn, 1, RealTime);
};

/**
 * Current based synapse with exponential decay.
 */
struct CurExp
    : public CurrentBasedSynapseBase<CurExp, CurExpState, CurExpParameters> {
	using Base = CurrentBasedSynapseBase<CurExp, CurExpState, CurExpParameters>;
	using Base::Base;
	using Base::p;

	friend SynapseBase<CurExp, CurExpState, CurExpParameters>;

private:
	Real m_tau_inv;

	template <typename State2, typename System>
	void process_spike(const Time &, Time, State2 &s, System &) const
	{
		s[0] += p().w_syn();
	}

public:
	/**
	 * Constructor of the CurExp synapse, allowing to directly set the synapse
	 * parameters.
	 *
	 * @param w_syn is the synaptic weight, determining by how much the
	 * current is increased whenever an external spike is received.
	 * @param tau_syn is the synaptic time constant, determining the slope of
	 * the exponential decay.
	 * @param input_spikes is a list containing the input spike trains.
	 */
	CurExp(Current w_syn, RealTime tau_syn,
	       const std::vector<Time> &input_spikes = std::vector<Time>())
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

	template <typename State2, typename System>
	State df(const State2 &s, const System &) const
	{
		// Exponential decay of the current
		return State({-m_tau_inv * s[0]});
	}
};
}

#endif /* CINDER_MODELS_SYNAPSES_CUR_EXP_HPP */

