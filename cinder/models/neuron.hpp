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
 * @file neuron.hpp
 *
 * Contains the base class which should be used for neuron implementations.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_NEURON_HPP
#define CINDER_MODELS_NEURON_HPP

#include <tuple>

#include <cinder/common/types.hpp>
#include <cinder/common/vector.hpp>
#include <cinder/ode/ode.hpp>

namespace cinder {

/**
 * Vector class representing both the neuron state and the current source state.
 */
template <typename MembraneState_, typename CurrentSourceState_>
class NeuronState
    : public VectorBase<NeuronState<MembraneState_, CurrentSourceState_>, Real,
                        MembraneState_::size() + CurrentSourceState_::size()> {
public:
	using MembraneState = MembraneState_;
	using CurrentSourceState = CurrentSourceState_;
	using Base =
	    VectorBase<NeuronState<MembraneState_, CurrentSourceState_>, Real,
	               MembraneState_::size() + CurrentSourceState_::size()>;
	using Base::Base;

	static constexpr NeuronState norm()
	{
		return concat<Real>(MembraneState::norm().as_array(),
		              CurrentSourceState::norm().as_array());
	}

	/**
	 * Returns a view at the neuron state vector.
	 */
	auto membrane_state()
	{
		return Base::template view<MembraneState::size(), 0>();
	}

	/**
	 * Returns a view at the neuron state vector.
	 */
	auto membrane_state() const
	{
		return Base::template view<MembraneState::size(), 0>();
	}

	/**
	 * Returns a view at the current source state vector.
	 */
	auto current_source_state()
	{
		return Base::template view<CurrentSourceState::size(),
		                           MembraneState::size()>();
	}

	/**
	 * Returns a view at the current source state vector.
	 */
	auto current_source_state() const
	{
		return Base::template view<CurrentSourceState::size(),
		                           MembraneState::size()>();
	}
};

/**
 * Implementation of a basic parameterised neuron membrane.
 */
template <typename State_, typename Parameters_, typename SpikeCallback_>
class MembraneBase : public ODEBase<State_> {
public:
	using Parameters = Parameters_;
	using SpikeCallback = SpikeCallback_;

private:
	Parameters m_params;

	SpikeCallback m_spike_callback;

protected:
	/**
	 * Member function to be called whenever the neuron implementation emits a
	 * spike.
	 */
	void emit_spike(Time t) { m_spike_callback(t); }

public:
	MembraneBase(const Parameters &params, const SpikeCallback &spike_callback)
	    : m_params(params), m_spike_callback(spike_callback)
	{
	}

	/**
	 * Returns the initial state vector. Sets the initial potential to the
	 * resting potential.
	 */
	State_ s0() const { return State_().v(p().v_rest()); }

	/**
	 * Returns a writable reference at the parameter vector.
	 */
	Parameters &p() { return m_params; }

	/**
	 * Returns a constant reference at the parameter vector.
	 */
	const Parameters &p() const { return m_params; }

	/**
	 * Retrieves the membrane potential from the current state vector.
	 */
	template <typename State, typename System>
	Voltage voltage(const State &s, const System &) const
	{
		return Voltage(s[0]);
	}
};

/**
 * Class which implements the spiking mechanism used in the LIF, AdEx and
 * Izhekevich models. Handles the generation and emission of spikes as well as
 * the refractory period.
 */
template <typename State_, typename Parameters_, typename SpikeCallback_,
          bool ArtificialSpike_>
class SpikingMembraneBase
    : public MembraneBase<State_, Parameters_, SpikeCallback_> {
private:
	using Base = MembraneBase<State_, Parameters_, SpikeCallback_>;
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
		const Voltage v_th = ArtificialSpike_ ? p().v_thresh() : p().v_spike();
		if (s[0] >= v_th) {
			// Plan the end of the refractory period
			m_ref_end = t + Time::sec(p().tau_refrac());

			if (ArtificialSpike_) {
				// Set the membrane potential to the spike potential
				s[0] = p().v_spike();
				sys.recorder().record(t, sys.s(), sys);
			}

			// Set the membrane potential to the reset potential
			s[0] = p().v_reset();
			emit_spike(t);
		}
	}
};

/**
 * Combines a membrane and a current source to a neuron.
 */
template <typename Membrane_, typename CurrentSource_>
class Neuron : public ODEBase<NeuronState<typename Membrane_::State,
                                          typename CurrentSource_::State>> {
public:
	using Membrane = Membrane_;
	using CurrentSource = CurrentSource_;
	using MembraneState = typename Membrane::State;
	using CurrentSourceState = typename CurrentSource::State;
	using State = NeuronState<MembraneState, CurrentSourceState>;

private:
	Membrane m_membrane;
	CurrentSource m_current_source;

public:
	Neuron(const Membrane &membrane, const CurrentSource &current_source)
	    : m_membrane(membrane), m_current_source(current_source)
	{
	}

	Membrane &membrane() { return m_membrane; }
	const Membrane &membrane() const { return m_membrane; }

	CurrentSource &current_source() { return m_current_source; }
	const CurrentSource &current_source() const { return m_current_source; }

	template <typename State, typename System>
	void init(Time t, const State &s, const System &sys)
	{
		m_membrane.init(t, s, sys);
		m_current_source.init(t, s, sys);
	}

	State s0() const
	{
		State res;
		res.membrane_state().assign(m_membrane.s0());
		res.current_source_state().assign(m_current_source.s0());
		return res;
	}

	Time next_discontinuity(Time t) const
	{
		return std::min(m_membrane.next_discontinuity(t),
		                m_current_source.next_discontinuity(t));
	}

	template <typename State, typename System>
	void handle_discontinuity(Time t, State &s, System &sys)
	{
		auto membrane_view = s.template view<MembraneState::size(), 0>();
		auto current_source_view = s.template view<CurrentSourceState::size(),
		                                           MembraneState::size()>();

		m_membrane.handle_discontinuity(t, membrane_view, sys);
		m_current_source.handle_discontinuity(t, current_source_view, sys);
	}

	template <typename State, typename System>
	void update(Time t, State &s, System &sys)
	{
		auto membrane_view = s.template view<MembraneState::size(), 0>();
		auto current_source_view = s.template view<CurrentSourceState::size(),
		                                           MembraneState::size()>();

		m_membrane.update(t, membrane_view, sys);
		m_current_source.update(t, current_source_view, sys);
	}

	template <typename State, typename System>
	Current current(const State &s, const System &sys) const
	{
		return m_current_source.current(
		    s.template view<CurrentSourceState::size(),
		                    MembraneState::size()>(),
		    sys);
	}

	template <typename State, typename System>
	Voltage voltage(const State &s, const System &sys) const
	{
		return m_membrane.voltage(s.template view<MembraneState::size(), 0>(),
		                          sys);
	}

	template <typename State2, typename System>
	State df(const State2 &s, const System &sys) const
	{
		State res;

		// Calculate the derivative of the neuron
		res.membrane_state().assign(
		    m_membrane.df(s.template view<MembraneState::size(), 0>(), sys));

		// Calculate the derivative of the current source
		res.current_source_state().assign(
		    m_current_source.df(s.template view<CurrentSourceState::size(),
		                                        MembraneState::size()>(),
		                        sys));

		return res;
	}
};

template <template <class> class Membrane, typename CurrentSource,
          typename SpikeCallback>
auto make_neuron(const CurrentSource &current_source,
                 const SpikeCallback &spike_callback,
                 const typename Membrane<SpikeCallback>::Parameters &params =
                     typename Membrane<SpikeCallback>::Parameters())
{
	return Neuron<Membrane<SpikeCallback>, CurrentSource>(
	    Membrane<SpikeCallback>(params, spike_callback), current_source);
}
}

#endif /* CINDER_MODELS_NEURON_HPP */
