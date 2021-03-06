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
 * @file izhikevich.hpp
 *
 * Implementation of the Izhikevich neuron model. Note that the currents
 * injected into the Izhikevich neuron under the assumption of a 1nF membrane
 * capacitance. In the original Izhikevich model the injected current has no
 * particular unit -- the capacitance is needed to convert between currents in
 * ampere and a voltage differential in volts per second.
 *
 * @author Andreas Stöckel
 * @author Christoph Jenzen
 */

#pragma once

#ifndef CINDER_MODELS_NEURONS_IZHIKEVICH_HPP
#define CINDER_MODELS_NEURONS_IZHIKEVICH_HPP

#include <cinder/models/neuron.hpp>

namespace cinder {
/**
 * The IzhikevichState class is the state vector which represents the state of
 * the Izhikevich membrane.
 */
struct IzhikevichState : public VectorBase<IzhikevichState, Real, 2> {
	using VectorBase<IzhikevichState, Real, 2>::VectorBase;

	TYPED_VECTOR_ELEMENT(v, 0, Voltage);
	NAMED_VECTOR_ELEMENT(u, 1);
};

/**
 * IzhikevichParameters is the parameter vector which contains the parameters
 * of the Izhikevich membrane. Note that the Izhikevich parameters are used
 * without any units. While considerably sloppy, this allows to directly use
 * parameter sets provided by Izhikevich and to stay compatible with PyNN.
 *
 * The membrane capcitance parameter cm is used to convert external currents to
 * volts per second.
 */
struct IzhikevichParameters : public VectorBase<IzhikevichParameters, Real, 5> {
	using VectorBase<IzhikevichParameters, Real, 5>::VectorBase;

	/**
	 * Default constructor. Initializes the neuron to the values also used in
	 * PyNN.
	 */
	IzhikevichParameters()
	{
		a(0.02);
		b(0.2);
		c(-65.0);
		d(2.0);
		cm(1_nF);
	}

	NAMED_VECTOR_ELEMENT(a, 0);
	NAMED_VECTOR_ELEMENT(b, 1);
	NAMED_VECTOR_ELEMENT(c, 2);
	NAMED_VECTOR_ELEMENT(d, 3);
	TYPED_VECTOR_ELEMENT(cm, 4, Capacitance);

	/**
	 * Returns the explicit refractory period -- the Izhikevich model does not
	 * possess an explicit refractory period, so zero is returned. Used by the
	 * SpikingMembraneBase class.
	 */
	static constexpr RealTime tau_refrac() { return RealTime(); }
	/**
	 * Returns the resting potential -- the Izhikevich model does not possess an
	 * explicit refractory period, so zero is returned. Used by the MembraneBase
	 * class.
	 */
	static constexpr Voltage v_rest() { return -70_mV; }
	/**
	 * Returns the reset potential in volt. Used by the SpikingMembraneBase
	 * class.
	 */
	Voltage v_reset() const { return Voltage(c() * 1e-3); }
	/**
	 * Returns the threshold potential in volt. The threshold potential is
	 * fixed in the Izhikevich model. Used by the SpikingMembraneBase class.
	 */
	static constexpr Voltage v_thresh() { return 30_mV; }
	/**
	 * Returns the spike potential in volt. The threshold potential is fixed in
	 * the Izhikevich model. Used by the SpikingMembraneBase class.
	 */
	static constexpr Voltage v_spike() { return v_thresh(); }
};

/**
 * Implementation fo the Izhikevich neuron model.
 */
class Izhikevich : public SpikingMembraneBase<Izhikevich, IzhikevichState,
                                              IzhikevichParameters, true> {
private:
	friend class SpikingMembraneBase<Izhikevich, IzhikevichState,
	                                 IzhikevichParameters, true>;

	using Base = SpikingMembraneBase<Izhikevich, IzhikevichState,
	                                 IzhikevichParameters, true>;

	Real m_cm_inv;

	template <typename State, typename System>
	void handle_output_spike(Time, State &s, System &)
	{
		s[1] += p().d();
	}

public:
	using Base::Base;
	using Base::p;

	template <typename State, typename System>
	void init(Time, const State &, const System &)
	{
		// Use the inverse of some values in order to avoid divisions in the
		// df code.
		m_cm_inv = 1.0_R / p().cm();
	}

	IzhikevichState s0() const
	{
		return IzhikevichState({p().v_rest(), -14.0_R});
	}

	template <typename State, typename System>
	IzhikevichState df(const State &s, const System &sys) const
	{
		const Real v = s[0] * 1e3_R;  // Convert V to mV
		const Real u = s[1];
		return IzhikevichState(
		    {(0.04_R * v * v + 5.0_R * v + 140.0_R - u) +
		         sys.ode().current(s, sys) * m_cm_inv,
		     p().a() * (p().b() * v - u) * 1e3_R});
	}
};
}

#endif /* CINDER_MODELS_NEURONS_IZHIKEVICH_HPP */	
