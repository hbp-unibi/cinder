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
 * Base class of all current sources. Provides the current() method.
 *
 * @tparam State_ is a vector which describes the current state of the current
 * source.
 * @tparam Parameters_ is a vector which contains the user-configurable
 * parameters of the current source.
 */
template <typename State_, typename Parameters_ = NullParameters>
struct CurrentSourceBase : public ODEBase<State_, Parameters_> {
	using ODEBase<State_, Parameters_>::ODEBase;

	/**
	 * Retrieves the current from the given state vector. Per default, the first
	 * state component is used. Current sources which do not possess any state
	 * (use the NullState) class must override this method, all other current
	 * source implementations are advised to place the current in the first
	 * state component (if possible).
	 */
	template <typename State, typename System>
	static Current current(const State &s, const System &)
	{
		return Current(s[0]);
	}
};

/**
 * Current source which injects no current. Can be used as a dummy current
 * source.
 */
struct NullCurrentSource : public CurrentSourceBase<NullState, NullParameters> {
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
 * ConstantCurrentSourceParameters contains the parameters describing a constant
 * current source.
 */
struct ConstantCurrentSourceParameters
    : public VectorBase<ConstantCurrentSourceParameters, Real, 1> {
	using VectorBase<ConstantCurrentSourceParameters, Real, 1>::VectorBase;

	TYPED_VECTOR_ELEMENT(i, 0, Current);
};

/**
 * Current source which injects a constant current into the neuron membrane.
 */
struct ConstantCurrentSource
    : public CurrentSourceBase<NullState, ConstantCurrentSourceParameters> {
	using Base = CurrentSourceBase<NullState, ConstantCurrentSourceParameters>;
	using Base::p;
	using Base::Base;

	ConstantCurrentSource(Current constant_current = 0.0_A)
	    : Base({{constant_current}})
	{
	}

	/**
	 * Returns the constant current stored in the parameters.
	 */
	template <typename State, typename System>
	Current current(const State &, const System &) const
	{
		return p().i();
	}
};

/**
 * State vector used by the StepCurrentSource class to describe the current that
 * is being injected into the neuron. The current i() is the current that is
 * currently being injected by the current source.
 */
struct StepCurrentSourceState
    : public VectorBase<StepCurrentSourceState, Real, 1> {
	using VectorBase<StepCurrentSourceState, Real, 1>::VectorBase;

	TYPED_VECTOR_ELEMENT(i, 0, Current);
};

/**
 * StepCurrentSourceParameters contains the parameters describing the behaviour
 * of the StepCurrentSource. There are three parameters: i(), which determines
 * the current that is being injected into the neuron, t_start() which
 * determines the time at which the current is first injected, and t_end() which
 * determines the time at which the current injection stops.
 */
struct StepCurrentSourceParameters
    : public VectorBase<StepCurrentSourceParameters, Real, 3> {
	using VectorBase<StepCurrentSourceParameters, Real, 3>::VectorBase;

	/**
	 * Constructor setting the default parameters.
	 */
	StepCurrentSourceParameters()
	{
		i(1_nA);
		t_start(0_s);
		t_end(RealTime(std::numeric_limits<Real>::max()));
	}

	TYPED_VECTOR_ELEMENT(i, 0, Current);
	TYPED_VECTOR_ELEMENT(t_start, 1, RealTime);
	TYPED_VECTOR_ELEMENT(t_end, 2, RealTime);
};

/**
 * A current source which injects a constant current during a given interval.
 */
class StepCurrentSource
    : public CurrentSourceBase<StepCurrentSourceState,
                               StepCurrentSourceParameters> {
private:
	Time m_t_start;
	Time m_t_end;

public:
	using Base =
	    CurrentSourceBase<StepCurrentSourceState, StepCurrentSourceParameters>;
	using Base::p;
	using Base::Base;

	StepCurrentSource(Current step_current, RealTime t_start,
	                  RealTime t_end = RealTime(std::numeric_limits<Real>::max()))
	    : Base({{step_current, t_start, t_end}})
	{
	}

	/**
	 * Converts the floating point times from the parameters to the internal
	 * fixed point times.
	 */
	template <typename State2, typename System>
	void init(Time, const State2 &, const System &)
	{
		m_t_start = Time::sec(p().t_start());
		m_t_end = Time::sec(p().t_end());
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
		s[0] = (t >= m_t_start && t < m_t_end) ? p().i() : 0_A;
	}
};

/**
 * A current source which may consist of a collection of current sources. This
 * allows to construct arbitrary synaptic input to the neuron. Use the
 * make_current_source function to conviniently construct a MultiCurrentSource
 * instance.
 */
template <typename... T>
class MultiCurrentSource : public MultiODE<T...> {
public:
	using Base = MultiODE<T...>;
	using Base::Base;

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

		return Base::template get<I>().current(
		           s.template view<InnerSize, Offs>(), sys) +
		       current_impl<State2, System, I + 1, Offs + InnerSize, Ts...>(
		           s, sys);
	}

public:
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
