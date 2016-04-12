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
 * @file vector.hpp
 *
 * Contains the implementation of a mathematical, fixed size vector suited for
 * automated parallelisation.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_COMMON_VECTOR_HPP
#define CINDER_COMMON_VECTOR_HPP

#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <type_traits>

namespace cinder {

/*
 * Forward declaration
 */
template <typename T, size_t Size>
class Vector;

/**
 * Base vector class without its own storage. Holds a reference to a segment of
 * an arbitrary array.
 */
template <typename Impl, typename T, typename Mem, size_t Size, size_t Offs>
class VectorView {
public:
	using own_type = VectorView<Impl, T, Mem, Size, Offs>;
	using value_type = T;
	using pointer = T *;
	using const_pointer = const T *;

private:
	/**
	 * Reference at the underlying storage.
	 */
	Mem &m_mem;

	template <typename Func, typename Other>
	friend Impl map(const own_type &v1, const Other &v2, Func f)
	{
		Impl res;
		for (size_t i = 0; i < Size; i++) {
			res[i] = f(v1[i], v2[i]);
		}
		return res;
	}

	template <typename Func>
	friend Impl map(const own_type &v, Func f)
	{
		Impl res;
		for (size_t i = 0; i < Size; i++) {
			res[i] = f(v[i]);
		}
		return res;
	}

	template <typename Func, typename Other>
	friend void map_onto(own_type &v1, const Other &v2, Func f)
	{
		for (size_t i = 0; i < Size; i++) {
			f(v1[i], v2[i]);
		}
	}

	template <typename Func>
	friend void map_onto(own_type &v, Func f)
	{
		for (size_t i = 0; i < Size; i++) {
			f(v[i]);
		}
	}

public:
	VectorView(Mem &mem) : m_mem(mem){};

	static constexpr size_t size() { return Size; }

	template <typename U = Mem, typename = typename std::enable_if<
	                                !std::is_const<U>::value>::type>
	pointer begin()
	{
		return &m_mem[Offs];
	}

	template <typename U = Mem, typename = typename std::enable_if<
	                                !std::is_const<U>::value>::type>
	pointer end()
	{
		return &m_mem[Offs + Size];
	}

	const_pointer begin() const { return &m_mem[Offs]; }

	const_pointer end() const { return &m_mem[Offs + Size]; }

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	void assign(const Other &v)
	{
		static_assert(own_type::size() == Other::size(),
		              "Vector size missmatch");
		std::copy(v.begin(), v.end(), begin());
	}

	void assign(value_type s)
	{
		map_onto(*this, [s](value_type &a) { return a = s; });
	}

	/**
	 * Returns a reference at the i-th element in the vector.
	 */
	value_type &operator[](size_t idx) { return m_mem[Offs + idx]; }

	/**
	 * Returns a copy of the i-th element in the vector.
	 */
	value_type operator[](size_t idx) const { return m_mem[Offs + idx]; }

	value_type sqrL2Norm() const
	{
		value_type res = 0;
		for (size_t i = 0; i < Size; i++) {
			res += (*this)[i] * (*this)[i];
		}
		return Size == 0 ? res : res * value_type(1.0 / double(Size));
	}

	value_type L2Norm() const { return std::sqrt(sqrL2Norm()); }

	friend std::ostream &operator<<(std::ostream &os, const own_type &m)
	{
		for (size_t i = 0; i < Size; i++) {
			os << (i == 0 ? "" : ", ") << m[i];
		}
		return os;
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend void operator+=(own_type &v1, const Other &v2)
	{
		static_assert(own_type::size() == Other::size(),
		              "Vector size missmatch");
		map_onto(v1, v2, [](value_type &a, value_type b) { return a += b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend void operator-=(own_type &v1, const Other &v2)
	{
		static_assert(own_type::size() == Other::size(),
		              "Vector size missmatch");
		map_onto(v1, v2, [](value_type &a, value_type b) { return a -= b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend void operator*=(own_type &v1, const Other &v2)
	{
		static_assert(own_type::size() == Other::size(),
		              "Vector size missmatch");
		map_onto(v1, v2, [](value_type &a, value_type b) { return a *= b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend void operator/=(own_type &v1, const Other &v2)
	{
		static_assert(own_type::size() == Other::size(),
		              "Vector size missmatch");
		map_onto(v1, v2, [](value_type &a, value_type b) { return a /= b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend Impl operator+(const own_type &v1, const Other &v2)
	{
		static_assert(own_type::size() == Other::size(),
		              "Vector size missmatch");
		return map(v1, v2, [](value_type a, value_type b) { return a + b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend Impl operator-(const own_type &v1, const Other &v2)
	{
		static_assert(own_type::size() == Other::size(),
		              "Vector size missmatch");
		return map(v1, v2, [](value_type a, value_type b) { return a - b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend Impl operator*(const own_type &v1, const Other &v2)
	{
		static_assert(own_type::size() == Other::size(),
		              "Vector size missmatch");
		return map(v1, v2, [](value_type a, value_type b) { return a * b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend Impl operator/(const own_type &v1, const Other &v2)
	{
		static_assert(own_type::size() == Other::size(),
		              "Vector size missmatch");
		return map(v1, v2, [](value_type a, value_type b) { return a / b; });
	}

	friend void operator+=(own_type &v, value_type s)
	{
		map_onto(v, [s](value_type &a) { return a += s; });
	}

	friend void operator-=(own_type &v, value_type s)
	{
		map_onto(v, [s](value_type &a) { return a -= s; });
	}

	friend void operator*=(own_type &v, value_type s)
	{
		map_onto(v, [s](value_type &a) { return a *= s; });
	}

	friend void operator/=(own_type &v, value_type s)
	{
		map_onto(v, [s](value_type &a) { return a /= s; });
	}

	friend Impl operator+(const own_type &v, value_type s)
	{
		return map(v, [s](value_type a) { return a + s; });
	}

	friend Impl operator-(const own_type &v, value_type s)
	{
		return map(v, [s](value_type a) { return a - s; });
	}

	friend Impl operator*(value_type s, const own_type &v)
	{
		return map(v, [s](value_type a) { return s * a; });
	}

	friend Impl operator*(const own_type &v, value_type s)
	{
		return map(v, [s](value_type a) { return a * s; });
	}

	friend Impl operator/(const own_type &v, value_type s)
	{
		return map(v, [s](value_type a) { return a / s; });
	}

	/**
	 * Returns a view onto this vector of ViewSize elements, starting at the
	 * element indicated by ViewOffs.
	 */
	template <size_t ViewSize, size_t ViewOffs = 0>
	auto view()
	{
		return VectorView<Vector<T, ViewSize>, T, Mem, ViewSize,
		                  Offs + ViewOffs>(m_mem);
	}

	template <size_t ViewSize, size_t ViewOffs = 0>
	auto view() const
	{
		return VectorView<Vector<T, ViewSize>, T, Mem, ViewSize,
		                  Offs + ViewOffs>(m_mem);
	}
};

template <typename Impl, typename T, size_t Size>
class alignas(16) VectorBase
    : public VectorView<Impl, T, std::array<T, Size>, Size, 0> {
public:
	using own_type = VectorBase<Impl, T, Size>;
	using base_type = VectorView<Impl, T, std::array<T, Size>, Size, 0>;
	using array_type = std::array<T, Size>;

private:
	array_type m_arr;

public:
	VectorBase() : VectorView<Impl, T, std::array<T, Size>, Size, 0>(m_arr) {}

	VectorBase(const own_type &o)
	    : VectorView<Impl, T, std::array<T, Size>, Size, 0>(m_arr),
	      m_arr(o.m_arr)
	{
		// Prevent VectorView from copying the memory pointer
	}

	own_type &operator=(const own_type &o)
	{
		// Prevent VectorView from copying the memory pointer
		m_arr = o.m_arr;
		return *this;
	}

	VectorBase(const array_type &arr)
	    : VectorView<Impl, T, std::array<T, Size>, Size, 0>(m_arr), m_arr(arr)
	{
	}

	VectorBase(std::initializer_list<T> init) : VectorBase()
	{
		std::copy(init.begin(), init.begin() + std::min(Size, init.size()),
		          m_arr.begin());
	}

	/**
	 * Returns a view onto this vector of ViewSize elements, starting at the
	 * element indicated by ViewOffs.
	 */
	template <size_t ViewSize, size_t ViewOffs = 0>
	auto view()
	{
		return VectorView<Vector<T, ViewSize>, T, std::array<T, Size>, ViewSize,
		                  ViewOffs>(m_arr);
	}

	template <size_t ViewSize, size_t ViewOffs = 0>
	auto view() const
	{
		return VectorView<Vector<T, ViewSize>, T, const std::array<T, Size>,
		                  ViewSize, ViewOffs>(m_arr);
	}
};

template <typename T, size_t Size>
class Vector final : public VectorBase<Vector<T, Size>, T, Size> {
public:
	using VectorBase<Vector, T, Size>::VectorBase;
};

#define NAMED_VECTOR_ELEMENT(NAME, IDX)           \
	static constexpr size_t idx_##NAME = IDX;     \
	void NAME(value_type x) { (*this)[IDX] = x; } \
	value_type &NAME() { return (*this)[IDX]; }   \
	value_type NAME() const { return (*this)[IDX]; }

#define TYPED_VECTOR_ELEMENT(NAME, IDX, TYPE)   \
	static constexpr size_t idx_##NAME = IDX;   \
	own_type& NAME(TYPE x) { (*this)[IDX] = x.v(); return *this;} \
	TYPE NAME() const { return TYPE((*this)[IDX]); }
}

#endif /* CINDER_COMMON_VECTOR_HPP */

