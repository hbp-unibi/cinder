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
 * @file mat.hpp
 *
 * Implementation of the multi-timescale adaptive threshold (MAT) model as
 * described in
 *
 * Made-to-order spiking neuron model equipped with a multi-timescale adaptive
 * threshold, Frontiers in Computational Neuroscience, 30 July 2009
 * Ryota Kobayashi, Yasuhiro Tsubo and Shigeru Shinomoto
 *
 * http://journal.frontiersin.org/article/10.3389/neuro.10.009.2009/full
 *
 * The MAT model can make use of an arbitrary number L of exponential functions
 * to model the threshold. Cinder uses an L + 1 dimensional state vector to
 * model these functions. This is required as it doesn't allow non-autonomous
 * differential equations.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_NEURONS_MAT_HPP
#define CINDER_MODELS_NEURONS_MAT_HPP

#include <cinder/models/neuron.hpp>

namespace cinder {

/**
 * The MATState class contains the current state of a MAT membrane with L = 2
 * independent time constants used for the exponentially decaying adaptive
 * threshold.
 */
struct MAT2State : public VectorBase<MAT2State, Real, 3> {
	using VectorBase<MAT2State, Real, 3>::VectorBase;

	TYPED_VECTOR_ELEMENT(v, 0, Voltage);
	TYPED_VECTOR_ELEMENT(th_1, 1, Voltage);
	TYPED_VECTOR_ELEMENT(th_2, 2, Voltage);
};

/**
 * The MAT2Parameters describes the parameters of a MAT neuron with L = 2
 * independent time constants describing the exponentially decaying adaptive
 * threshold.
 */
struct MAT2Parameters : public VectorBase<MAT2Parameters, Real, 10> {
	using VectorBase<MAT2Parameters, Real, 10>::VectorBase;

	/**
	 * Default constructor.
	 */
	MAT2Parameters()
	{
		cm(0.5_nF);
		tau_m(5_ms);
		v_thresh(-46_mV);
		v_rest(-65_mV);
		v_spike(20_mV);
		tau_refrac(2_ms);
		alpha_1(37_ms);
		alpha_2(2_ms);
		tau_1(10_ms);
		tau_2(200_ms);
	}

	TYPED_VECTOR_ELEMENT(cm, 0, Capacitance);
	TYPED_VECTOR_ELEMENT(g_leak, 1, Conductance);
	TYPED_VECTOR_ELEMENT(v_thresh, 2, Voltage);
	TYPED_VECTOR_ELEMENT(v_rest, 3, Voltage);
	TYPED_VECTOR_ELEMENT(v_spike, 4, Voltage);
	TYPED_VECTOR_ELEMENT(tau_refrac, 5, RealTime);
	TYPED_VECTOR_ELEMENT(alpha_1, 6, RealTime);
	TYPED_VECTOR_ELEMENT(alpha_2, 7, RealTime);
	TYPED_VECTOR_ELEMENT(tau_1, 8, RealTime);
	TYPED_VECTOR_ELEMENT(tau_2, 9, RealTime);

	MAT2Parameters &tau_m(RealTime t)
	{
		g_leak(Conductance(cm().v() / t.v()));
		return *this;
	}

	RealTime tau_m() const { return RealTime(cm().v() / g_leak().v()); }
};

/**
 * Class implementing the multi-timescale adaptive threshold model described in
 *
 * Made-to-order spiking neuron model equipped with a multi-timescale adaptive
 * threshold, Frontiers in Computational Neuroscience, 30 July 2009
 * Ryota Kobayashi, Yasuhiro Tsubo and Shigeru Shinomoto
 *
 * http://journal.frontiersin.org/article/10.3389/neuro.10.009.2009/full
 */
class MAT2 : public MembraneBase<MAT2State, MAT2Parameters> {
private:
	using Base = MembraneBase<MAT2State, MAT2Parameters>;

	Real m_cm_inv;
	Real m_tau_1_inv;
	Real m_tau_2_inv;
	Time m_ref_end = MAX_TIME;

public:
	using Base::Base;
	using Base::p;
	using Base::emit_spike;

	/**
	 * Returns true if the membrane currently is in its refractory state.
	 */
	bool in_refrac() const { return m_ref_end != MAX_TIME; }

	/**
	 * Returns the time at which the refractory period will end.
	 */
	Time next_discontinuity(Time) const { return m_ref_end; }

	/**
	 * Ends the refractory period once its end is reached.
	 */
	template <typename State, typename System>
	void handle_discontinuity(Time t, State &, System &)
	{
		if (m_ref_end != MAX_TIME && t > m_ref_end) {
			m_ref_end = MAX_TIME;
		}
	}

	/**
	 * Called just before the simulation is started. Calculates the inverse of
	 * the membrane capacitance in order to keep the costly division out of the
	 * df() method.
	 */
	template <typename State, typename System>
	void init(Time, const State &, const System &)
	{
		m_cm_inv = 1.0_R / p().cm();
		m_tau_1_inv = 1.0_R / p().tau_1();
		m_tau_2_inv = 1.0_R / p().tau_2();
	}

	/**
	 * Calculates the voltage differential for the MAT membrane given the
	 * current state and an external current.
	 *
	 * @param s is the current neuron state.
	 * @param sys is a reference at the global system state, allowing to access
	 * the current that is being injected into the neuron.
	 * @return a LIFState vector containing the time-differential of the state.
	 */
	template <typename State, typename System>
	MAT2State df(const State &s, const System &sys) const
	{
		const Current i_syn{sys.ode().current(s, sys)};
		const Current i_rest{(p().v_rest() - s[0]) * p().g_leak()};
		return MAT2State({
			(i_rest + i_syn) * m_cm_inv,
			-s[1] * m_tau_1_inv,
			-s[2] * m_tau_2_inv
		});
	}

	/**
	 * Emit a spike whenever the dynamic threshold potential is surpassed.
	 */
	template <typename State, typename System>
	void update(Time t, State &s, System &sys)
	{
		// Check whether the current membrane potential crosses the threshold
		// potential
		const Real v = s[0];
		const Real v_th = p().v_thresh() + s[1] + s[2];
		if (v > v_th && !in_refrac()) {
			// Plan the end of the refractory period
			const RealTime tau_refrac = p().tau_refrac();
			if (tau_refrac > 0) {
				m_ref_end = t + Time::sec(tau_refrac);
			}

			// Increase the threshold by increasing the individual threshold
			// components
			s[1] += p().alpha_1();
			s[2] += p().alpha_2();

			// Temporarily set the membrane potential to the spike potential
			// and record the current system state
			s[0] = p().v_spike();
			sys.recorder().record(t, sys.s(), sys, true);
			s[0] = v;

			// Emit the spike
			emit_spike(t);
		}
	}
};
}

#endif /* CINDER_MODELS_NEURONS_MAT_HPP */
