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
 * @file time.hpp
 *
 * Contains the definition of the Time type. To avoid non-uniform resolution in
 * Time deltas -- as it would occur with floating point values for time -- the
 * Time type is implemented as a 64-bit integer with fixed point arithmetic.
 * Time can be converted to the RealTime type from types.hpp, which is a
 * floating point time type.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_COMMON_TIME_HPP
#define CINDER_COMMON_TIME_HPP

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>

namespace cinder {

/**
 * Integer type used internally by the Time type to represent times.
 */
using TimeType = int64_t;

/**
 * Factor for converting a floating point time in seconds to a time value.
 */
constexpr double SEC_TO_TIME = double(1L << 48);

/**
 * Factor for converting a time value to seconds.
 */
constexpr double TIME_TO_SEC = 1.0 / SEC_TO_TIME;

/**
 * Maximum internal time value.
 */
constexpr TimeType MAX_INT_TIME = std::numeric_limits<TimeType>::max();

/**
 * Minimum internal time value.
 */
constexpr TimeType MIN_INT_TIME = std::numeric_limits<TimeType>::min();

/**
 * Time is the type used for storing time values. Times are stored as a 64 bit
 * integer in with fixed divisor 2^48. This is done to avoid artifacts, e.g.
 * when calculating time deltas, where the resolution would decrease with
 * increasing time for floating point values. This header also defines a set
 * of suffixes, which allow to create time values with the corresponding SI
 * suffixes.
 *
 * @seealso RealTime a floating-point version of Time. Time can be implicitly
 * converted to RealTime, but not vice-versa.
 */
struct Time {
private:
	/**
	 * Static version of the "fromSeconds" method.
	 *
	 * @param ft is the floating point value for which the corresponding
	 * integer value should be returned.
	 * @return the internal integer value used to represent the time.
	 */
	static constexpr TimeType secondsToTimeType(double ft)
	{
		return (int64_t)(ft * SEC_TO_TIME) > MAX_INT_TIME
		           ? MAX_INT_TIME
		           : ((int64_t)(ft * SEC_TO_TIME) < MIN_INT_TIME
		                  ? MIN_INT_TIME
		                  : ft * SEC_TO_TIME);
	}

public:
	/**
	 * Internal fixed-point integer Time value.
	 */
	TimeType t;

	/**
	 * Default constructor, initializes the time value to zero.
	 */
	Time() : t(0) {}

	/**
	 * Constructor, initializes the tiem value with the given internal integer
	 * timestamp.
	 *
	 * @param t is the internal integer time the Time instance should be
	 * initialized to.
	 * @seealso Time::sec creates a new Time instance from a floating point
	 * value which denotes a time in seconds.
	 * @seealso Time::msec creates a new Time instance from a floating point
	 * value which denotes a time in milliseconds.
	 */
	explicit constexpr Time(TimeType t) : t(t) {}

	/**
	 * Creates a new instance of the Time structure from a float.
	 *
	 * @param t is a time in seconds.
	 */
	static constexpr Time sec(double t) { return Time(secondsToTimeType(t)); }

	/**
	 * Creates a new instance of the Time structure from a float.
	 *
	 * @param t is a time in milliseconds.
	 */
	static constexpr Time msec(double t)
	{
		return Time(secondsToTimeType(t * 1e-3));
	}

	/**
	 * Converts the internal integer Time to a floating point time in seconds.
	 *
	 * @return the current time value in seconds.
	 */
	constexpr double sec() const { return double(t) * TIME_TO_SEC; }

	/**
	 * Converts the internal integer Time to a floating point time in
	 * milliseconds.
	 *
	 * @return the current time value in milliseconds.
	 */
	constexpr double msec() const { return double(t) * TIME_TO_SEC * 1e3; }

	/* Operators */

	friend Time operator+(const Time &t1) { return Time(t1.t); }
	friend Time operator-(const Time &t1) { return Time(-t1.t); }
	friend Time operator+(const Time &t1, const Time &t2)
	{
		return Time(t1.t + t2.t);
	}
	friend Time operator-(const Time &t1, const Time &t2)
	{
		return Time(t1.t - t2.t);
	}
	friend Time operator/(const Time &t, double s)
	{
		return Time(t.t / s);
	}
	friend Time operator*(const Time &t, double s)
	{
		return Time(TimeType(t.t * s));
	}
	friend Time operator*(double s, const Time &t)
	{
		return Time(TimeType(s * t.t));
	}
	friend Time operator%(const Time &t1, const Time &t2)
	{
		return Time(t1.t % t2.t);
	}
	friend void operator+=(Time &t1, const Time &t2) { t1.t += t2.t; }
	friend void operator-=(Time &t1, const Time &t2) { t1.t -= t2.t; }
	friend void operator*=(Time &t, double s) { t.t *= s; }
	friend void operator/=(Time &t, double s) { t.t /= s; }
	friend bool operator==(const Time &t1, const Time &t2)
	{
		return t1.t == t2.t;
	}
	friend bool operator!=(const Time &t1, const Time &t2)
	{
		return t1.t != t2.t;
	}
	friend bool operator<(const Time &t1, const Time &t2)
	{
		return t1.t < t2.t;
	}
	friend bool operator<=(const Time &t1, const Time &t2)
	{
		return t1.t <= t2.t;
	}
	friend bool operator>(const Time &t1, const Time &t2)
	{
		return t1.t > t2.t;
	}
	friend bool operator>=(const Time &t1, const Time &t2)
	{
		return t1.t >= t2.t;
	}

	/**
	 * Prints the time in seconds to the given output stream.
	 *
	 * @param os is the outuput stream to which the time value should be written.
	 * @param t is the time value that should be printed.
	 */
	friend std::ostream &operator<<(std::ostream &os, const Time &t)
	{
		return os << t.sec();
	}
};

static constexpr Time operator"" _s(long double t)
{
	return Time::sec(t * 1e0);
}
static constexpr Time operator"" _s(unsigned long long int t)
{
	return Time::sec(t * 1e0);
}
static constexpr Time operator"" _ms(long double t)
{
	return Time::sec(t * 1e-3);
}
static constexpr Time operator"" _ms(unsigned long long int t)
{
	return Time::sec(t * 1e-3);
}
static constexpr Time operator"" _us(long double t)
{
	return Time::sec(t * 1e-6);
}
static constexpr Time operator"" _us(unsigned long long int t)
{
	return Time::sec(t * 1e-6);
}
static constexpr Time operator"" _ns(long double t)
{
	return Time::sec(t * 1e-9);
}
static constexpr Time operator"" _ns(unsigned long long int t)
{
	return Time::sec(t * 1e-9);
}
static constexpr Time operator"" _ps(long double t)
{
	return Time::sec(t * 1e-12);
}
static constexpr Time operator"" _ps(unsigned long long int t)
{
	return Time::sec(t * 1e-12);
}
static constexpr Time operator"" _fs(long double t)
{
	return Time::sec(t * 1e-15);
}
static constexpr Time operator"" _fs(unsigned long long int t)
{
	return Time::sec(t * 1e-15);
}

/**
 * Maximum representable time.
 */
constexpr Time MAX_TIME = Time(MAX_INT_TIME);

/**
 * Maximum representable time.
 */
constexpr Time MIN_TIME = Time(MIN_INT_TIME);

/**
 * Maximum representable time in seconds.
 */
constexpr double MAX_TIME_SEC = MAX_INT_TIME / SEC_TO_TIME;

/**
 * Maximum representable time in seconds.
 */
constexpr double MIN_TIME_SEC = MIN_INT_TIME / SEC_TO_TIME;

/**
 * Minimum time difference representable by the Time type.
 */
constexpr double MIN_TIME_DELTA = TIME_TO_SEC;

}

namespace std {
/**
 * Overload for allowing to apply std::abs to the Time type.
 */
inline cinder::Time abs(cinder::Time t) { return (t.t < 0) ? -t : t; }
}

#endif /* CINDER_COMMON_TIME_HPP */

