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
 * @file hodgkin_huxley.hpp
 *
 * Implementation of the four dimensional Hodgkin-Huxley neuron model.
 *
 * @author Andreas Stöckel
 * @author Christoph Jenzen
 */

#pragma once

#ifndef CINDER_MODELS_NEURONS_HODGKIN_HUXLEY_HPP
#define CINDER_MODELS_NEURONS_HODGKIN_HUXLEY_HPP

#include <cmath>

#include <cinder/models/neuron.hpp>

namespace cinder {
/**
 * The HodgkinHuxleyState class is the state vector which represents the state
 * of the HodgkinHuxley membrane.
 */
struct HodgkinHuxleyState : public VectorBase<HodgkinHuxleyState, Real, 4> {
	using VectorBase<HodgkinHuxleyState, Real, 4>::VectorBase;

	TYPED_VECTOR_ELEMENT(v, 0, Voltage);
	NAMED_VECTOR_ELEMENT(n, 1);
	NAMED_VECTOR_ELEMENT(m, 2);
	NAMED_VECTOR_ELEMENT(h, 3);
};

/**
 * HodgkinHuxleyParameters is the parameter vector which contains the parameters
 * of the HodgkinHuxley membrane.
 */
struct HodgkinHuxleyParameters
    : public VectorBase<HodgkinHuxleyParameters, Real, 8> {
	using VectorBase<HodgkinHuxleyParameters, Real, 8>::VectorBase;

	/**
	 * Default constructor. Initializes the neuron to the values also used in
	 * PyNN and the hh_cond_exp_traub class.
	 */
	HodgkinHuxleyParameters()
	{
		gbar_Na(20_uS);
		gbar_K(6_uS);
		g_leak(0.01_uS);
		cm(0.2_nF);
		v_offset(-63_mV);
		e_rev_Na(50_mV);
		e_rev_K(-90_mV);
		v_rest(-65_mV);
	}

	TYPED_VECTOR_ELEMENT(gbar_Na, 0, Conductance);
	TYPED_VECTOR_ELEMENT(gbar_K, 1, Conductance);
	TYPED_VECTOR_ELEMENT(g_leak, 2, Conductance);
	TYPED_VECTOR_ELEMENT(cm, 3, Capacitance);
	TYPED_VECTOR_ELEMENT(v_offset, 4, Voltage);
	TYPED_VECTOR_ELEMENT(e_rev_leak, 5, Voltage);
	TYPED_VECTOR_ELEMENT(e_rev_Na, 6, Voltage);
	TYPED_VECTOR_ELEMENT(e_rev_K, 7, Voltage);

	HodgkinHuxleyParameters &tau_m(RealTime t)
	{
		g_leak(Conductance(cm().v() / t.v()));
		return *this;
	}

	RealTime tau_m() const { return RealTime(cm().v() / g_leak().v()); }

	/**
	 * Alias for the e_rev_leak parameter.
	 */
	Voltage v_rest() const { return e_rev_leak(); }

	/**
	 * Returns the current resting potential, alias for the e_rev_leak
	 * parameter.
	 */
	HodgkinHuxleyParameters &v_rest(Voltage v)
	{
		e_rev_leak(v);
		return *this;
	}

	/**
	 * Returns the explicit refractory period -- the HodgkinHuxley model does
	 * not possess an explicit refractory period, so zero is returnex. Used by
	 * the SpikingMembraneBase class.
	 */
	static constexpr RealTime tau_refrac() { return RealTime(); }
};

/**
 * Implementation fo the Hodgkin-Huxley neuron model with generic channel
 * dynamics.
 */
template <typename ChannelDynamics>
class HodgkinHuxleyBase
    : public MembraneBase<HodgkinHuxleyState, HodgkinHuxleyParameters> {
private:
	using Base = MembraneBase<HodgkinHuxleyState, HodgkinHuxleyParameters>;

	/**
	 * The EvaluatedChannelDynamics helper class is used to evaluated the alpha
	 * and beta functions provided by the ChannelDynamics type.
	 */
	struct EvaluatedChannelDynamics {
		const Real alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h;

		EvaluatedChannelDynamics(const Real v)
		    : alpha_n(ChannelDynamics::alpha_n(v)),
		      beta_n(ChannelDynamics::beta_n(v)),
		      alpha_m(ChannelDynamics::alpha_m(v)),
		      beta_m(ChannelDynamics::beta_m(v)),
		      alpha_h(ChannelDynamics::alpha_h(v)),
		      beta_h(ChannelDynamics::beta_h(v))
		{
		}
	};

	Voltage m_last_v;
	bool m_in_refrac;

	static Real pow2(Real x) { return x * x; }
	static Real pow3(Real x) { return pow2(x) * x; }
	static Real pow4(Real x) { return pow2(pow2(x)); }

public:
	using Base::Base;
	using Base::p;
	using Base::emit_spike;

	template <typename State, typename System>
	void init(Time, const State &s, const System &)
	{
		m_in_refrac = false;
		m_last_v = Voltage(s[0]) - p().v_offset();
	}

	HodgkinHuxleyState s0() const
	{
		EvaluatedChannelDynamics x((p().e_rev_leak() - p().v_offset()) * 1e3_R);
		return HodgkinHuxleyState({p().v_rest(),
		                           x.alpha_n / (x.alpha_n + x.beta_n),
		                           x.alpha_m / (x.alpha_m + x.beta_m),
		                           x.alpha_h / (x.alpha_h + x.beta_h)});
	}

	template <typename State, typename System>
	HodgkinHuxleyState df(const State &s, const System &sys) const
	{
		const Real n = s[HodgkinHuxleyState::idx_n];
		const Real m = s[HodgkinHuxleyState::idx_m];
		const Real h = s[HodgkinHuxleyState::idx_h];

		const Real n4 = pow4(s[HodgkinHuxleyState::idx_n]);
		const Real m3 = pow3(s[HodgkinHuxleyState::idx_m]);

		const Real i_Na = m3 * h * p().gbar_Na() * (s[0] - p().e_rev_Na());
		const Real i_K = n4 * p().gbar_K() * (s[0] - p().e_rev_K());
		const Real i_l = p().g_leak() * (s[0] - p().e_rev_leak());
		const Real i_syn = sys.ode().current(s, sys);
		const Real dtV = (i_syn - (i_Na + i_K + i_l)) / p().cm();

		const Real v = (s[0] - p().v_offset()) * 1e3_R;  // Convert V to mV
		const EvaluatedChannelDynamics x(v);
		const Real dtN = (x.alpha_n - (x.alpha_n + x.beta_n) * n) * 1e3_R;
		const Real dtM = (x.alpha_m - (x.alpha_m + x.beta_m) * m) * 1e3_R;
		const Real dtH = (x.alpha_h - (x.alpha_h + x.beta_h) * h) * 1e3_R;

		return HodgkinHuxleyState({dtV, dtN, dtM, dtH});
	}

	template <typename State, typename System>
	void update(Time t, const State &s, const System &)
	{
		const Voltage v = Voltage(s[0]) - p().v_offset();
		if (v >= 30_mV) {
			if (!m_in_refrac && v < m_last_v) {
				m_in_refrac = true;
				emit_spike(t);
			}
		}
		else {
			m_in_refrac = false;
		}
		m_last_v = v;
	}
};

/**
 * Structure describing the Hodgkin-Huxley channel dynamics as described in
 *
 * Traub, R.x. and Miles, R. (1991)
 * Neuronal Networks of the Hippocampus. Cambridge University Press,
 * Cambridge UK.
 *
 * The formulas were directly taken from the NEST source-code of the
 * hh_cond_exp_traub model.
 */
class TraubChannelDynamics {
private:
	static Real limit_denom(Real x)
	{
		static constexpr Real MIN_DENOM = 1e-6;
		if (std::abs(x) < MIN_DENOM) {
			if (x < 0) {
				return -MIN_DENOM;
			}
			else {
				return MIN_DENOM;
			}
		}
		return x;
	}

	static Real limit(Real x)
	{
		static constexpr Real MAX_RES = 1e3;
		return std::min(MAX_RES, std::max(-MAX_RES, x));
	}

public:
	static Real alpha_n(Real V)
	{
		return limit(0.032_R * (15.0_R - V) /
		             limit_denom(std::exp((15._R - V) / 5._R) - 1._R));
	}

	static Real beta_n(Real V)
	{
		return limit(0.5_R * std::exp((10._R - V) / 40._R));
	}

	static Real alpha_m(Real V)
	{
		return limit(0.32_R * (13._R - V) /
		             limit_denom((std::exp((13._R - V) / 4._R) - 1._R)));
	}

	static Real beta_m(Real V)
	{
		return limit(0.28_R * (V - 40._R) /
		             limit_denom((std::exp((V - 40._R) / 5._R) - 1._R)));
	}

	static Real alpha_h(Real V)
	{
		return limit(0.128_R * std::exp((17._R - V) / 18._R));
	}

	static Real beta_h(Real V)
	{
		return limit(4._R / limit_denom((1._R + std::exp((40._R - V) / 5._R))));
	}
};

/**
 * Implementation of a Hodgkin-Huxley type neuron with TraubChannelDynamics.
 */
class HodgkinHuxley : public HodgkinHuxleyBase<TraubChannelDynamics> {
public:
	using HodgkinHuxleyBase<TraubChannelDynamics>::HodgkinHuxleyBase;
};
}

#endif /* CINDER_MODELS_NEURONS_HODGKIN_HUXLEY_HPP */
