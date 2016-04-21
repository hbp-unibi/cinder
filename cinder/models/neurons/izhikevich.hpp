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
 * Implementation of the Izhikevich neuron model.
 *
 * @author Andreas Stöckel
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

	static constexpr IzhikevichState norm()
	{
		return IzhikevichState({1e3, 1e0});
	}
};

/**
 * IzhikevichParameters is the parameter vector which contains the parameters
 * of the Izhikevich membrane. Note that the Izhikevich parameters are used 
 * without any units. While considerably sloppy, this allows to directly use
 * parameter sets provided by Izhikevich and to stay compatible with PyNN.
 */
struct IzhikevichParameters : public VectorBase<IzhikevichParameters, Real, 4> {
	using VectorBase<IzhikevichParameters, Real, 4>::VectorBase;

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
	}

	NAMED_VECTOR_ELEMENT(a, 0);
	NAMED_VECTOR_ELEMENT(b, 1);
	NAMED_VECTOR_ELEMENT(c, 2);
	NAMED_VECTOR_ELEMENT(d, 3);

	/**
	 * Returns the explicit refractory period -- the Izhikevich model does not
	 * possess an explicit refractory period, so zero is returned. Used by the
	 * SpikingMembraneBase class.
	 */
	static constexpr RealTime tau_refrac() { return RealTime(); }

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
};

/**
 * Implementation fo the Izhikevich neuron model.
 */
template <typename SpikeCallback_>
class Izhikevich
    : public SpikingMembraneBase<Izhikevich<SpikeCallback_>, IzhikevichState,
                                 IzhikevichParameters, SpikeCallback_, false> {
private:
	using Base =
	    SpikingMembraneBase<Izhikevich<SpikeCallback_>, IzhikevichState,
	                        IzhikevichParameters, SpikeCallback_, false>;

	template<typename State, typename System>
	void handle_output_spike(Time, State &s, System &) {
		s[1] += p().d();
	}

public:
	using Base::Base;
	using Base::p;

	template <typename State, typename System>
	IzhikevichState df(const State &s, const System &sys) const
	{
		const Real v = s[0] * Real(1e3);  // Convert V to mV
		const Real u = s[1];
		return IzhikevichState(
		    {(Real(0.04) * v * v + Real(5.0) * v + Real(140.0) - u) *
		             Real(1e-3) +
		         sys.ode().current(s, sys),
		     p().a() * (p().b() * v - u)});
	}
};
}

#endif /* CINDER_MODELS_NEURONS_IZHIKEVICH_HPP */
