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

struct LIFState : public VectorBase<LIFState, Real, 1> {
	using VectorBase<LIFState, Real, 1>::VectorBase;
};

struct LIFParameters : public VectorBase<LIFParameters, Real, 7> {
	using VectorBase<LIFParameters, Real, 7>::VectorBase;

	LIFParameters()
	{
		cM(1_nF);
		gL(2_uS);
		eTh(-54_mV);
		eRest(-70_mV);
		eReset(-80_mV);
		eSpike(20_mV);
		tauRef(1_ms);
	}

	TYPED_VECTOR_ELEMENT(cM, 0, Capacitance);
	TYPED_VECTOR_ELEMENT(gL, 1, Conductance);
	TYPED_VECTOR_ELEMENT(eTh, 2, Voltage);
	TYPED_VECTOR_ELEMENT(eRest, 3, Voltage);
	TYPED_VECTOR_ELEMENT(eReset, 4, Voltage);
	TYPED_VECTOR_ELEMENT(eSpike, 5, Voltage);
	TYPED_VECTOR_ELEMENT(tauRef, 6, RealTime);
};

template <typename SpikeCallback_>
class LIF: public MembraneBase<LIFState, LIFParameters, SpikeCallback_> {
private:
	using Base = MembraneBase<LIFState, LIFParameters, SpikeCallback_>;
	Time m_ref_end = MAX_TIME;

public:
	using Base::Base;
	using Base::p;
	using Base::emit_spike;

	LIFState s0() const {
		return LIFState({p().eRest()});
	}

	/**
	 * Returns the time at which the refractory period will end.
	 */
	Time next_discontinuity(Time t) const
	{
		return t >= m_ref_end ? MAX_TIME : m_ref_end;
	}


	template <typename State, typename System>
	void handle_discontinuity(Time t, State &, System &) {
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
		if (s[0] > p().eTh()) {
			// Plan the end of the refractory period
			m_ref_end = t + Time::sec(p().tauRef());

			// Set the membrane potential to the spike potential
			s[0] = p().eSpike();
			sys.recorder().record(t, sys.s(), sys);

			// Set the membrane potential to the reset potential
			s[0] = p().eReset();
			emit_spike(t);
		}
	}

	template <typename State, typename System>
	LIFState df(const State &s, const System &sys) const
	{
		if (m_ref_end == MAX_TIME) {
			const Current i_syn{sys.ode().current(s, sys)};
			const Current i_rest{(p().eRest() - s[0]) * p().gL()};
			return LIFState({(i_rest + i_syn) / p().cM()});
		} else {
			return LIFState({0});
		}
	}
};
}

#endif /* CINDER_MODELS_NEURONS_LIF_HPP */
