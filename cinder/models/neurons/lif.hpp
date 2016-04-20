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
 * @file lif.hpp
 *
 * Implementation of the simple linear integrate and fire neuron model.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_NEURONS_LIF_HPP
#define CINDER_MODELS_NEURONS_LIF_HPP

#include <cinder/models/neuron.hpp>

namespace cinder {
/**
 * The LIFState class is the state vector which represents the state of the
 * LIF membrane.
 */
struct LIFState : public VectorBase<LIFState, Real, 1> {
	using VectorBase<LIFState, Real, 1>::VectorBase;

	static constexpr LIFState norm() { return LIFState({1e3}); }
};

/**
 * LIFParameters is the parameter vector which contains the state of the LIF
 * membrane.
 */
struct LIFParameters : public VectorBase<LIFParameters, Real, 7> {
	using VectorBase<LIFParameters, Real, 7>::VectorBase;

	/**
	 * Default constructor. Initializes the neuron to the values also used in
	 * PyNN.
	 */
	LIFParameters()
	{
		cm(1_nF);
		tau_m(20_ms);
		v_thresh(-50_mV);
		v_rest(-65_mV);
		v_reset(-65_mV);
		v_spike(20_mV);
		tau_refrac(0.1_ms);
	}

	TYPED_VECTOR_ELEMENT(cm, 0, Capacitance);
	TYPED_VECTOR_ELEMENT(g_leak, 1, Conductance);
	TYPED_VECTOR_ELEMENT(v_thresh, 2, Voltage);
	TYPED_VECTOR_ELEMENT(v_rest, 3, Voltage);
	TYPED_VECTOR_ELEMENT(v_reset, 4, Voltage);
	TYPED_VECTOR_ELEMENT(v_spike, 5, Voltage);
	TYPED_VECTOR_ELEMENT(tau_refrac, 6, RealTime);

	LIFParameters &tau_m(RealTime t)
	{
		g_leak(Conductance(cm().v() / t.v()));
		return *this;
	}

	RealTime tau_m() const { return RealTime(cm().v() / g_leak().v()); }
};

template <typename SpikeCallback_>
class LIF : public MembraneBase<LIFState, LIFParameters, SpikeCallback_> {
private:
	using Base = MembraneBase<LIFState, LIFParameters, SpikeCallback_>;
	Time m_ref_end = MAX_TIME;

public:
	using Base::Base;
	using Base::p;
	using Base::emit_spike;

	LIFState s0() const { return LIFState({p().v_rest()}); }

	/**
	 * Returns the time at which the refractory period will end.
	 */
	Time next_discontinuity(Time) const { return m_ref_end; }

	template <typename State, typename System>
	void handle_discontinuity(Time t, State &, System &)
	{
		if (t >= m_ref_end) {
			m_ref_end = MAX_TIME;
		}
	}

	/**
	 * Emit a spike whenever the threshold potential is surpassed.
	 */
	template <typename State, typename System>
	void update(Time t, State &s, System &sys)
	{
		if (s[0] > p().v_thresh()) {
			// Plan the end of the refractory period
			m_ref_end = t + Time::sec(p().tau_refrac());

			// Set the membrane potential to the spike potential
			s[0] = p().v_spike();
			sys.recorder().record(t, sys.s(), sys);

			// Set the membrane potential to the reset potential
			s[0] = p().v_reset();
			emit_spike(t);
		}
	}

	template <typename State, typename System>
	LIFState df(const State &s, const System &sys) const
	{
		if (m_ref_end == MAX_TIME) {
			const Current i_syn{sys.ode().current(s, sys)};
			const Current i_rest{(p().v_rest() - s[0]) * p().g_leak()};
			return LIFState({(i_rest + i_syn) / p().cm()});
		}
		else {
			return LIFState({0});
		}
	}
};
}

#endif /* CINDER_MODELS_NEURONS_LIF_HPP */
