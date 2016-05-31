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
 * Implementation of the simple linear integrate-and-fire neuron model, see for
 * example chapter one of
 *
 * Neuronal Dynamics: From Single Neurons to Networks and Models of Cognition
 * Wulfram Gerstner, Werner M. Kistler, Richard Naud, Liam Paninski
 * Cambridge University Press, 2014
 *
 * http://neuronaldynamics.epfl.ch/online/Ch1.S3.html
 *
 * or alternatively chapter eight of
 *
 * Dynamical Systems in Neuroscience: The Geometry of Excitability and Bursting
 * Eugene M. Izhikevich
 * The MIT Press, 2007
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
 * LIF membrane. Contains a single element, namely the membrane voltage of the
 * LIF membrane.
 */
struct LIFState : public VectorBase<LIFState, Real, 1> {
	using VectorBase<LIFState, Real, 1>::VectorBase;

	TYPED_VECTOR_ELEMENT(v, 0, Voltage);
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

/**
 * Class implementing the classical linear integrate-and-fire model. This class
 * is derived from the "SpikingMembraneBase" class which implements the
 * production of output spikes, the reset behaviour and refractoriness.
 */
class LIF : public SpikingMembraneBase<LIF, LIFState, LIFParameters, true> {
private:
	using Base = SpikingMembraneBase<LIF, LIFState, LIFParameters, true>;

	/**
	 * Inverse of the membrane capacitance -- calculated in the init function in
	 * order to keep the costly division out of the df() method.
	 */
	Real m_cm_inv;

public:
	using Base::Base;
	using Base::in_refrac;
	using Base::p;

	/**
	 * Called just before the simulation is started. Calculates the inverse of
	 * the membrane capacitance in order to keep the costly division out of the
	 * df() method.
	 */
	template <typename State, typename System>
	void init(Time, const State &, const System &)
	{
		m_cm_inv = 1.0_R / p().cm();
	}

	/**
	 * Calculates the voltage differential for the LIF membrane given the
	 * current state and an external current.
	 *
	 * @param s is the current neuron state.
	 * @param sys is a reference at the global system state, allowing to access
	 * the current that is being injected into the neuron.
	 * @return a LIFState vector containing the time-differential of the state.
	 */
	template <typename State, typename System>
	LIFState df(const State &s, const System &sys) const
	{
		if (in_refrac()) {
			return LIFState();  // If in refractory period, return a zero-vector
		}
		const Current i_syn{sys.ode().current(s, sys)};
		const Current i_rest{(p().v_rest() - s[0]) * p().g_leak()};
		return LIFState({(i_rest + i_syn) * m_cm_inv});
	}
};
}

#endif /* CINDER_MODELS_NEURONS_LIF_HPP */
