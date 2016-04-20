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
 * @file ode.hpp
 *
 * Contains a basic implementation of the ODE concept.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_ODE_ODE_HPP
#define CINDER_ODE_ODE_HPP

#include <tuple>

#include <cinder/common/time.hpp>
#include <cinder/common/types.hpp>
#include <cinder/common/vector.hpp>

namespace cinder {
/**
 * The ODEBase represents a basic ordinary differential equation with no
 * discontinuities.
 *
 * @tparam State is the state vector underlying the ODE.
 */
template <typename State_>
struct ODEBase {
	using State = State_;

	/**
	 * Allows the implementations to perform some initialization.
	 */
	template <typename State, typename System>
	static void init(Time, const State &, const System &)
	{
	}

	/**
	 * Returns the initial state of the current source.
	 */
	static State s0() { return State(); }

	/**
	 * Returns the position of the next discontinuity in the differential
	 * equation.
	 */
	static Time next_discontinuity(Time) { return MAX_TIME; }

	/**
	 * Called whenever the discontinuity is reached.
	 */
	template <typename State, typename System>
	void handle_discontinuity(Time, const State &, const System &)
	{
	}

	/**
	 * Gives the model the opportunity to perform some additional calculations.
	 */
	template <typename State, typename System>
	static void update(Time, const State &, const System &)
	{
	}

	/**
	 * Returns the differential.
	 */
	template <typename State, typename System>
	static State_ df(const State &, const System &)
	{
		return State_();
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

/* See
 * http://stackoverflow.com/questions/25068481/c11-constexpr-flatten-list-of-stdarray-into-array
 */

template <size_t... Is>
struct seq {
};
template <size_t N, size_t... Is>
struct gen_seq : gen_seq<N - 1, N - 1, Is...> {
};
template <size_t... Is>
struct gen_seq<0, Is...> : seq<Is...> {
};

template <typename T, size_t N1, size_t... I1, size_t N2, size_t... I2>
// Expansion pack
constexpr std::array<T, N1 + N2> concat(const std::array<T, N1> &a1,
                                        const std::array<T, N2> &a2, seq<I1...>,
                                        seq<I2...>)
{
	return {a1[I1]..., a2[I2]...};
}

template <typename T, size_t N1, size_t N2>
// Initializer for the recursion
constexpr std::array<T, N1 + N2> concat(const std::array<T, N1> &a1,
                                        const std::array<T, N2> &a2)
{
	return concat(a1, a2, gen_seq<N1>{}, gen_seq<N2>{});
}
}

/**
 * Class allowing to concatenate multiple ODE classes into one.
 */
template <typename... T>
class MultiODE {
public:
	/**
	 * State vector representing the entire system.
	 */
	class State : public VectorBase<State, Real, internal::dim<T...>()> {
	private:
		template <size_t I>
		static constexpr auto norm_impl()
		{
			return std::array<Real, 0>();
		}

		template <size_t I, typename T0, typename... Ts>
		static constexpr auto norm_impl()
		{
			return internal::concat<Real>(T0::State::norm().as_array(),
			                              norm_impl<0, Ts...>());
		}

	public:
		using VectorBase<State, Real, internal::dim<T...>()>::VectorBase;

		static constexpr State norm()
		{
			return State(norm_impl<0, T...>());
		}
	};

protected:
	std::tuple<T...> instances;

private:
	/*
	 * Compile-time recursive implementation of init.
	 */

	template <typename State2, typename System, size_t I, size_t Offs>
	void init_impl(Time, const State2 &, const System &) const
	{
	}

	template <typename State2, typename System, size_t I, size_t Offs,
	          typename T0, typename... Ts>
	void init_impl(Time t, const State2 &s, const System &sys) const
	{
		using InnerState = typename T0::State;
		static constexpr size_t InnerSize = InnerState::size();

		get<I>().init(t, s.template view<InnerSize, Offs>(), sys);
		init_impl<State2, System, I + 1, Offs + InnerSize, Ts...>(t, s, sys);
	}

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

		InnerState tmp = get<I>().s0();
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
		return std::min(get<I>().next_discontinuity(t),
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
		get<I>().handle_discontinuity(t, s_view, sys);

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

		auto state_view = s.template view<InnerSize, Offs>();
		get<I>().update(t, state_view, sys);

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

		auto tmp = get<I>().df(s.template view<InnerSize, Offs>(), sys);
		std::copy(tmp.begin(), tmp.end(), res.begin() + Offs);

		df_impl<State2, System, I + 1, Offs + InnerSize, Ts...>(res, s, sys);
	}

public:
	MultiODE(const T &... args) : instances(args...) {}

	template <typename State2, typename System>
	void init(Time t, const State2 &s, const System &sys)
	{
		init_impl<State2, System, 0, 0, T...>(t, s, sys);
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

	/**
	 * Returns the initial state of the current source.
	 */
	State s0() const
	{
		State res;
		s0_impl<0, 0, T...>(res);
		return res;
	}

	/**
	 * Returns the position of the next discontinuity in the differential
	 * equation.
	 */
	Time next_discontinuity(Time t) const
	{
		return next_discontinuity_impl<0, T...>(t);
	}

	/**
	 * Called whenever the discontinuity is reached.
	 */
	template <typename State2, typename System>
	void handle_discontinuity(Time t, State2 &s, System &sys)
	{
		handle_discontinuity_impl<State2, System, 0, 0, T...>(t, s, sys);
	}

	/**
	 * Returns the differential.
	 */
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
};
}

#endif /* CINDER_ODE_ODE_HPP */
