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
 * @file adex.hpp
 *
 * Implementation of the adaptive exponential integrate and fire neuron.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_NEURONS_ADEX_HPP
#define CINDER_MODELS_NEURONS_ADEX_HPP

#include <cinder/common/fast_math.hpp>
#include <cinder/models/neuron.hpp>

namespace cinder {
/**
 * The AdExState class is the state vector which represents the state of the
 * AdExState membrane.
 */
struct AdExState : public VectorBase<AdExState, Real, 2> {
	using VectorBase<AdExState, Real, 2>::VectorBase;

	TYPED_VECTOR_ELEMENT(v, 0, Voltage);
	TYPED_VECTOR_ELEMENT(w, 1, Current);

	static constexpr AdExState scale() { return AdExState({1e3, 1e9}); }
};

/**
 * LIFParameters is the parameter vector which contains the state of the LIF
 * membrane.
 */
struct AdExParameters : public VectorBase<AdExParameters, Real, 11> {
	using VectorBase<AdExParameters, Real, 11>::VectorBase;

	/**
	 * Default constructor. Initializes the neuron to the values also used in
	 * PyNN.
	 */
	AdExParameters()
	{
		cm(0.281_nF);
		tau_m(9.3667_ms);
		v_thresh(-50.4_mV);
		v_rest(-70.6_mV);
		v_reset(-70.6_mV);
		v_spike(-40.0_mV);
		tau_refrac(0.1_ms);
		delta_T(2.0_mV);
		a(4.0_nS);
		b(0.0805_nA);
		tau_w(144.0_ms);
	}

	TYPED_VECTOR_ELEMENT(cm, 0, Capacitance);
	TYPED_VECTOR_ELEMENT(g_leak, 1, Conductance);
	TYPED_VECTOR_ELEMENT(v_thresh, 2, Voltage);
	TYPED_VECTOR_ELEMENT(v_rest, 3, Voltage);
	TYPED_VECTOR_ELEMENT(v_reset, 4, Voltage);
	TYPED_VECTOR_ELEMENT(v_spike, 5, Voltage);
	TYPED_VECTOR_ELEMENT(tau_refrac, 6, RealTime);
	TYPED_VECTOR_ELEMENT(delta_T, 7, Voltage);
	TYPED_VECTOR_ELEMENT(a, 8, Conductance);
	TYPED_VECTOR_ELEMENT(b, 9, Current);
	TYPED_VECTOR_ELEMENT(tau_w, 10, RealTime);

	AdExParameters &tau_m(RealTime t)
	{
		g_leak(Conductance(cm().v() / t.v()));
		return *this;
	}

	RealTime tau_m() const { return RealTime(cm().v() / g_leak().v()); }
};

template <typename SpikeCallback_>
class AdEx : public SpikingMembraneBase<AdEx<SpikeCallback_>, AdExState,
                                        AdExParameters, SpikeCallback_, false> {
private:
	friend class SpikingMembraneBase<AdEx<SpikeCallback_>, AdExState,
	                                 AdExParameters, SpikeCallback_, false>;
	using Base = SpikingMembraneBase<AdEx<SpikeCallback_>, AdExState,
	                                 AdExParameters, SpikeCallback_, false>;
	Real m_max_i_th_exp;
	Real m_delta_t_inv;
	Real m_cm_inv;
	Real m_tau_m_inv;
	Real m_tau_w_inv;

	/**
	 * Function called whenever an output spike is generated. Allows some models
	 * to update their interal state whenever an output spike occurs.
	 */
	template<typename State, typename System>
	void handle_output_spike(Time, State &s, System &) {
		s[1] += p().b();
	}

public:
	using Base::Base;
	using Base::in_refrac;
	using Base::p;

	template <typename State, typename System>
	void init(Time, const State &, const System &)
	{
		// Calculate the maximum exponent allowed in the threshold current
		static constexpr RealTime MIN_DELTA_T = 0.1_us;
		m_max_i_th_exp = log((p().v_spike() - p().v_reset()) /
		                   (MIN_DELTA_T.v() * p().delta_T() / p().tau_m()));

		// Use the inverse of some values in order to avoid divisions in the
		// df code.
		m_delta_t_inv = 1.0 / p().delta_T();
		m_cm_inv = 1.0 / p().cm();
		m_tau_w_inv = 1.0 / p().tau_w();
		m_tau_m_inv = 1.0 / p().tau_m();
	}

	template <typename State, typename System>
	AdExState df(const State &s, const System &sys) const
	{
		const Real dw = (p().a() * (s[0] - p().v_rest()) - s[1]) * m_tau_w_inv;
		if (in_refrac()) {
			return AdExState({0, dw});
		}
		const Current i_syn{sys.ode().current(s, sys)};
		const Current i_rest{(p().v_rest() - s[0]) * p().g_leak()};
		const Current i_th_exp{(s[0] - p().v_thresh()) * m_delta_t_inv};
		const Current i_th{p().g_leak().v() * p().delta_T().v() *
		                   fast::exp(std::min(m_max_i_th_exp, i_th_exp.v()))};
		const Current i_w{s[1]};
		const Real dv = (i_rest + i_syn + i_th - i_w) * m_cm_inv;
		return AdExState({dv, dw});
	}
};
}

#endif /* CINDER_MODELS_NEURONS_ADEX_HPP */
