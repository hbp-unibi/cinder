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

#include <tuple>

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
};

/**
 * State vector used by current sources with a single current component.
 */
struct SingleCurrentState : public VectorBase<SingleCurrentState, Real, 1> {
	using VectorBase<SingleCurrentState, Real, 1>::VectorBase;
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

namespace internal {
/**
 * Class used to calculate the total dimensionality of a state vector.
 */
template <typename... T>
struct Dim;

template <>
struct Dim<> {
	static constexpr size_t calculate() { return 0; }
};

template <typename H, typename... T>
struct Dim<H, T...> {
	static constexpr size_t calculate()
	{
		using State = typename H::State;
		return State::size() + Dim<T...>::calculate();
	}
};

template <typename... T>
static constexpr size_t dim()
{
	return Dim<T...>::calculate();
}
}

/**
 * A current source which may consist of a collection of current sources. This
 * allows to construct arbitrary synaptic input to the neuron. Use the
 * make_current_source function to conviniently construct a MultiCurrentSource
 * instance.
 */
template <typename... T>
class MultiCurrentSource {
public:
	/**
	 * State vector representing the entire system.
	 */
	class State : public VectorBase<State, Real, internal::dim<T...>()> {
		using VectorBase<State, Real, internal::dim<T...>()>::VectorBase;
	};

private:
	std::tuple<T...> instances;

	/*
	 * Compile-time recursive implementation of s0.
	 */
	template <size_t I, size_t Offs>
	void s0_impl(State &) const
	{
	}

	template <size_t I, size_t Offs, typename T0, typename... Ts>
	void s0_impl(State &res) const
	{
		using InnerState = typename T0::State;
		static constexpr size_t InnerSize = InnerState::size();

		InnerState tmp = std::get<I>(instances).s0();
		std::copy(tmp.begin(), tmp.begin() + InnerSize, res.begin() + Offs);

		s0_impl<I + 1, Offs + InnerSize, Ts...>(res);
	}

	/*
	 * Compile-time recursive implementation of next_discontinuity.
	 */

	template <size_t I>
	Time next_discontinuity_impl(Time) const
	{
		return MAX_TIME;
	}

	template <size_t I, typename T0, typename... Ts>
	Time next_discontinuity_impl(Time t) const
	{
		return std::min(std::get<I>(instances).next_discontinuity(t),
		                next_discontinuity_impl<I + 1, Ts...>(t));
	}

	/*
	 * Compile-time recursive implementation of handle_discontinuity.
	 */

	template <typename State2, typename System, size_t I, size_t Offs>
	void handle_discontinuity_impl(Time, State2 &, System &)
	{
	}

	template <typename State2, typename System, size_t I, size_t Offs,
	          typename T0, typename... Ts>
	void handle_discontinuity_impl(Time t, State2 &s, System &sys)
	{
		using InnerState = typename T0::State;
		static constexpr size_t InnerSize = InnerState::size();

		auto s_view = s.template view<InnerSize, Offs>();
		std::get<I>(instances).handle_discontinuity(t, s_view, sys);

		handle_discontinuity_impl<State2, System, I + 1, Offs + InnerSize,
		                          Ts...>(t, s, sys);
	}

	/*
	 * Compile-time recursive implementation of update.
	 */

	template <typename State2, typename System, size_t I, size_t Offs>
	void update_impl(Time, State2 &, System &)
	{
	}

	template <typename State2, typename System, size_t I, size_t Offs,
	          typename T0, typename... Ts>
	void update_impl(Time t, State2 &s, System &sys)
	{
		using InnerState = typename T0::State;
		static constexpr size_t InnerSize = InnerState::size();

		std::get<I>(instances)
		    .update(t, s.template view<InnerSize, Offs>(), sys);

		update_impl<State2, System, I + 1, Offs + InnerSize, Ts...>(t, s, sys);
	}

	/*
	 * Compile-time recursive implementation of df.
	 */

	template <typename State2, typename System, size_t I, size_t Offs>
	void df_impl(State &, const State2 &, const System &) const
	{
	}

	template <typename State2, typename System, size_t I, size_t Offs,
	          typename T0, typename... Ts>
	void df_impl(State &res, const State2 &s, const System &sys) const
	{
		using InnerState = typename T0::State;
		static constexpr size_t InnerSize = InnerState::size();

		auto tmp =
		    std::get<I>(instances).df(s.template view<InnerSize, Offs>(), sys);
		std::copy(tmp.begin(), tmp.end(), res.begin() + Offs);

		df_impl<State2, System, I + 1, Offs + InnerSize, Ts...>(res, s, sys);
	}

	/*
	 * Compile-time recursive implementation of current.
	 */
	template <typename State2, typename System, size_t I, size_t Offs>
	Current current_impl(const State2 &, const System &) const
	{
		return Current(0.0);
	}

	template <typename State2, typename System, size_t I, size_t Offs, typename T0,
	          typename... Ts>
	Current current_impl(const State2 &s, const System &sys) const
	{
		using InnerState = typename T0::State;
		static constexpr size_t InnerSize = InnerState::size();

		return std::get<I>(instances).current(s.template view<InnerSize, Offs>(), sys) +
		       current_impl<State2, System, I + 1, Offs + InnerSize, Ts...>(s, sys);
	}

public:
	MultiCurrentSource(const T &... args) : instances(args...) {}

	State s0() const
	{
		State res;
		s0_impl<0, 0, T...>(res);
		return res;
	}

	template <size_t I>
	auto &get()
	{
		return std::get<I>(instances);
	}

	template <size_t I>
	const auto &get() const
	{
		return std::get<I>(instances);
	}

	Time next_discontinuity(Time t) const
	{
		return next_discontinuity_impl<0, T...>(t);
	}

	template <typename State2, typename System>
	void handle_discontinuity(Time t, State2 &s, System &sys)
	{
		handle_discontinuity_impl<State2, System, 0, 0, T...>(t, s, sys);
	}

	template <typename State2, typename System>
	State df(const State2 &s, const System &sys) const
	{
		State res;
		df_impl<State2, System, 0, 0, T...>(res, s, sys);
		return res;
	}

	/**
	 * Gives the model the opportunity to perform some additional calculations.
	 */
	template <typename State2, typename System>
	void update(Time t, State2 &s, System &sys)
	{
		update_impl<State2, System, 0, 0, T...>(t, s, sys);
	}

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
