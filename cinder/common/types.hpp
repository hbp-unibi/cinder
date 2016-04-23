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
 * @file types.hpp
 *
 * Contains type definitions used throughout Cinder.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_COMMON_TYPES_HPP
#define CINDER_COMMON_TYPES_HPP

#include <cinder/config.h>
#include <cinder/common/time.hpp>

namespace cinder {

#if !defined(CINDER_REAL_WIDTH)
#define CINDER_REAL_WIDTH 8
#endif

/**
 * The Real type is a floating point type. Its width can be chosen by the user
 * by setting the CINDER_REAL_WIDTH macro.
 */
#if CINDER_REAL_WIDTH == 4
using Real = float;
#elif CINDER_REAL_WIDTH == 8
using Real = double;
#elif CINDER_REAL_WIDTH == 10
using Real = long double;
#elif CINDER_REAL_WIDTH == 16
extern "C" {
#include <quadmath.h>
}
using Real = __float128;
#else
using Real = float;
#error Invalid value for CINDER_REAL_WIDTH supplied!
#endif

/**
 * Type to be used in suffix declarations.
 */
using sx_double_t = long double;

/**
 * Type to be used in suffix declarations.
 */
using sx_int_t = unsigned long long int;

/**
 * The quantity type is a class wrapper for any C++ numerical type. This allows
 * to specialise numerical types, which is useful when encoding units in the
 * type system.
 */
template <typename Impl, typename T>
struct Quantity {
private:
	using Self = Quantity<Impl, T>;

	/**
	 * The actual value stored in this instance.
	 */
	T m_val;

public:
	/**
	 * Default constructor. Default-initializes the encapsulated value.
	 */
	constexpr Quantity() : m_val(T()) {}

	/**
	 * Creates a new instance of Quantity for the given value.
	 */
	explicit constexpr Quantity(T val) : m_val(val) {}

	/**
	 * Allows to implicitly convert Quantity to the underlying value type.
	 */
	constexpr operator T() const { return m_val; }

	constexpr T v() const { return m_val; }

	friend constexpr Time operator+(const Self &q) { return Time(q.v()); }
	friend constexpr Impl operator-(const Self &q) { return Impl(-q.v()); }

	friend constexpr Impl operator+(const Self &q1, const Self &q2)
	{
		return Impl(q1.v() + q2.v());
	}

	constexpr Impl& operator+=(const Self &q) {
		m_val += q.v();
		return static_cast<Impl&>(*this);
	}

	constexpr Impl& operator-=(const Self &q) {
		m_val -= q.v();
		return static_cast<Impl&>(*this);
	}

	friend constexpr Impl operator-(const Self &q1, const Self &q2)
	{
		return Impl(q1.v() - q2.v());
	}

	friend constexpr T operator/(const Self &q1, const Self &q2)
	{
		return T(q1.v() / q2.v());
	}

	friend constexpr Impl operator/(const Self &q, T s)
	{
		return Impl(q.v() / s);
	}
	friend constexpr Impl operator*(const Self &q, T s)
	{
		return Impl(q.v() * s);
	}
	friend constexpr Impl operator*(T s, const Self &q)
	{
		return Impl(s * q.v());
	}

	bool constexpr operator==(const Self &o) const { return v() == o.v(); }
	bool constexpr operator!=(const Self &o) const { return v() != o.v(); }
	bool constexpr operator<(const Self &o) const { return v() < o.v(); }
	bool constexpr operator<=(const Self &o) const { return v() <= o.v(); }
	bool constexpr operator>(const Self &o) const { return v() > o.v(); }
	bool constexpr operator>=(const Self &o) const { return v() >= o.v(); }
};

/**
 * Unit suffixes for currents
 */
struct Current : public Quantity<Current, Real> {
	using Quantity<Current, Real>::Quantity;
};

struct Voltage : public Quantity<Voltage, Real> {
	using Quantity<Voltage, Real>::Quantity;
};

struct Conductance : public Quantity<Conductance, Real> {
	using Quantity<Conductance, Real>::Quantity;
};

struct Capacitance : public Quantity<Capacitance, Real> {
	using Quantity<Capacitance, Real>::Quantity;
};

struct RealTime : public Quantity<RealTime, Real> {
	using Quantity<RealTime, Real>::Quantity;
	constexpr RealTime() : Quantity<RealTime, Real>() {}
	constexpr RealTime(Time t) : Quantity<RealTime, Real>(t.sec()) {}
};

#define DEFINE_SUFFIX(CLASS, PREFIX, SUFFIX, FAC)                      \
	static constexpr CLASS operator"" _##PREFIX##SUFFIX(sx_double_t x) \
	{                                                                  \
		return CLASS(x * FAC);                                         \
	}                                                                  \
	static constexpr CLASS operator"" _##PREFIX##SUFFIX(sx_int_t x)    \
	{                                                                  \
		return CLASS(x * FAC);                                         \
	}

#define DEFINE_SUFFIX_SET(CLASS, PREFIX)   \
	DEFINE_SUFFIX(CLASS, f, PREFIX, 1e-15) \
	DEFINE_SUFFIX(CLASS, p, PREFIX, 1e-12) \
	DEFINE_SUFFIX(CLASS, n, PREFIX, 1e-9)  \
	DEFINE_SUFFIX(CLASS, u, PREFIX, 1e-6)  \
	DEFINE_SUFFIX(CLASS, m, PREFIX, 1e-3)  \
	DEFINE_SUFFIX(CLASS, c, PREFIX, 1e-2)  \
	DEFINE_SUFFIX(CLASS, d, PREFIX, 1e-1)  \
	DEFINE_SUFFIX(CLASS, , PREFIX, 1e0)    \
	DEFINE_SUFFIX(CLASS, da, PREFIX, 1e1)  \
	DEFINE_SUFFIX(CLASS, h, PREFIX, 1e2)   \
	DEFINE_SUFFIX(CLASS, k, PREFIX, 1e3)   \
	DEFINE_SUFFIX(CLASS, M, PREFIX, 1e6)   \
	DEFINE_SUFFIX(CLASS, G, PREFIX, 1e9)   \
	DEFINE_SUFFIX(CLASS, T, PREFIX, 1e12)  \
	DEFINE_SUFFIX(CLASS, P, PREFIX, 1e15)

DEFINE_SUFFIX_SET(Current, A)
DEFINE_SUFFIX_SET(Voltage, V)
DEFINE_SUFFIX_SET(Conductance, S)
DEFINE_SUFFIX_SET(Capacitance, F)
}

#endif /* CINDER_COMMON_TYPES_HPP */

