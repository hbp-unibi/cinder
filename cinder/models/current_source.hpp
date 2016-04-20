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
 * @file current_source.hpp
 *
 * Implementation of some basic current sources. Furthermore implements the
 * MultiCurrentSource object which allows to combine an arbitrary number of
 * current sources into a single current source.
 *
 * A current source in itself is an ODE object which can be integrated using
 * the Solver class.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_MODELS_CURRENT_SOURCE_HPP
#define CINDER_MODELS_CURRENT_SOURCE_HPP

#include <cinder/common/types.hpp>
#include <cinder/common/time.hpp>
#include <cinder/common/vector.hpp>
#include <cinder/ode/ode.hpp>

namespace cinder {
/**
 * A current source which injects a zero current into the neuron.
 *
 * @tparam State is the neuron
 */
template <typename StateImpl>
struct CurrentSourceBase : public ODEBase<StateImpl> {
	/**
	 * Retrieves the current from the given state vector.
	 */
	template <typename State, typename System>
	static Current current(const State &s, const System &)
	{
		return Current(s[0]);
	}
};

/**
 * State vector used by current sources with no state component.
 */
struct NullState : public VectorBase<NullState, Real, 0> {
	using VectorBase<NullState, Real, 0>::VectorBase;

	static constexpr NullState norm() { return NullState(); }
};

/**
 * State vector used by current sources with a single current component.
 */
struct SingleCurrentState : public VectorBase<SingleCurrentState, Real, 1> {
	using VectorBase<SingleCurrentState, Real, 1>::VectorBase;

	static constexpr SingleCurrentState norm()
	{
		return SingleCurrentState({1e9});
	}
};

/**
 * Current source which injects zero current.
 */
struct NullCurrentSource : public CurrentSourceBase<NullState> {
	/**
	 * Retrieves the current from the given state vector.
	 */
	template <typename State, typename System>
	static Current current(const State &, const System &)
	{
		return 0_A;
	}
};

/**
 * Current source which injects a constant current.
 */
struct ConstantCurrentSource : public CurrentSourceBase<SingleCurrentState> {
private:
	Current m_constant_current;

public:
	ConstantCurrentSource(Current constant_current = 0.0_A)
	    : m_constant_current(constant_current)
	{
	}

	SingleCurrentState s0() const
	{
		return SingleCurrentState({m_constant_current});
	}
};

/**
 * A current source which injects a constant current during a given interval.
 */
class StepCurrentSource : public CurrentSourceBase<SingleCurrentState> {
private:
	Current m_step_current;
	Time m_t_start;
	Time m_t_end;

public:
	StepCurrentSource(Current m_step_current, Time m_t_start,
	                  Time m_t_end = MAX_TIME)
	    : m_step_current(m_step_current), m_t_start(m_t_start), m_t_end(m_t_end)
	{
	}

	/**
	 * Returns the position of the next discontinuity in the differential
	 * equation.
	 */
	Time next_discontinuity(Time t) const
	{
		return (t < m_t_start ? m_t_start : (t < m_t_end ? m_t_end : MAX_TIME));
	}

	/**
	 * Called whenever the discontinuity is reached.
	 */
	template <typename State, typename System>
	void handle_discontinuity(Time t, State &s, System &)
	{
		s[0] = (t >= m_t_start && t < m_t_end) ? m_step_current : 0_A;
	}
};

/**
 * A current source which may consist of a collection of current sources. This
 * allows to construct arbitrary synaptic input to the neuron. Use the
 * make_current_source function to conviniently construct a MultiCurrentSource
 * instance.
 */
template <typename... T>
class MultiCurrentSource: public MultiODE<T...> {
public:
	using Base = MultiODE<T...>;

private:
	/*
	 * Compile-time recursive implementation of current.
	 */
	template <typename State2, typename System, size_t I, size_t Offs>
	Current current_impl(const State2 &, const System &) const
	{
		return Current(0.0);
	}

	template <typename State2, typename System, size_t I, size_t Offs,
	          typename T0, typename... Ts>
	Current current_impl(const State2 &s, const System &sys) const
	{
		using InnerState = typename T0::State;
		static constexpr size_t InnerSize = InnerState::size();

		return Base::template get<I>()
		           .current(s.template view<InnerSize, Offs>(), sys) +
		       current_impl<State2, System, I + 1, Offs + InnerSize, Ts...>(
		           s, sys);
	}

public:
	using Base::Base;

	/**
	 * Retrieves the current from the given state vector.
	 */
	template <typename State2, typename System>
	Current current(const State2 &s, const System &sys) const
	{
		return current_impl<State2, System, 0, 0, T...>(s, sys);
	}
};

/**
 * Convenience function which creates a new MultiCurrentSource instance.
 */
template <typename... T>
MultiCurrentSource<T...> make_current_source(const T &... args)
{
	return MultiCurrentSource<T...>(args...);
}
}

#endif /* CINDER_MODELS_CURRENT_SOURCE_HPP */
