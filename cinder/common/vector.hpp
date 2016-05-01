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

#include <cinder/common/array_utils.hpp>

namespace cinder {

/*
 * Forward declarations
 */
template <typename, size_t>
class Vector;

template <typename T, size_t>
class VectorView;

/**
 * The VectorElementInfo structure stores information about the entry of a
 * vector. For example, the Vector class is used to represent the neuron
 * parameters -- the VectorElementInfo struct is then used to access the unit
 * and name of each parameter at runtime.
 */
struct VectorElementInfo {
	/**
	 * Name of the vector element. In a parameter vector a name could for
	 * example be "v_rest" for the resting potential of a neuron membrane.
	 */
	const char *name;

	/**
	 * Physical unit of the vector element. Should be an empty string if this
	 * is a unit-less element. Value should be without the SI prefix. For
	 * example, for a value describing a voltage, set unit to "V".
	 */
	const char *unit;

	/**
	 * Typical scale factor the values should be multiplied with to have the
	 * values in a unit-range without requiring an exponential. For example, if
	 * the values are mostly in a millivolt range, set scale to 1e3.
	 */
	double scale;

	/**
	 * Constructor of the VectorElementInfo class.
	 *
	 * @param name is the name of the vector element.
	 * @param unit is the physical unit of the vector element (if present)
	 * without SI suffix.
	 * @param scale is a scale factors the values should be multiplied with
	 * to get to a human-readable range without requiring exponential notation.
	 */
	constexpr VectorElementInfo(const char *name = "", const char *unit = "",
	                            double scale = 1e0)
	    : name(name), unit(unit), scale(scale)
	{
	}
};

/**
 * The VectorMixin class defines basic functionality that should be provided by
 * a constant-sized mathematical vector, such as addition, subtraction,
 * iteration, assignment operations and calculation of the L2 norm. This class
 * does not allocate any memory for the storage of the vector elements, it relys
 * on an implementation class Impl to provide a pointer at the memory.
 *
 * @tparam Impl is the type of the class which derives from VectorMixin and
 * which provides access at the underlying storage.
 * @tparam Inst is the type of vector instance that should be created as the
 * result of an operator.
 * @tparam T is the underlying storage type.
 * @tparam Size_ is the number of elements stored in this vector.
 */
template <typename Impl, typename Inst, typename T, size_t Size_>
class VectorMixin {
public:
	/**
	 * Constant containing the size of the vector.
	 */
	static constexpr size_t Size = Size_;

	/**
	 * Type alias for the complete type of the VectorMixin.
	 */
	using Self = VectorMixin<Impl, Inst, T, Size_>;

	/**
	 * Type of the element stored in the vector.
	 */
	using value_type = T;

	/**
	 * Pointer type of the element stored in the vector.
	 */
	using pointer = T *;

	/**
	 * Const pointer type of the element stored in the vector.
	 */
	using const_pointer = const T *;

private:
	/**
	 * Function which uses the Impl class to get access at the underlying memory
	 * pointer.
	 *
	 * @return a pointer at the first element represented by this vector.
	 */
	T *mem() { return static_cast<Impl &>(*this).mem(); }

	/**
	 * Function which uses the Impl class to get access at the underlying const
	 * memory pointer.
	 *
	 * @return a pointer at the first element represented by this vector.
	 */
	const T *mem() const { return static_cast<const Impl &>(*this).mem(); }

	/**
	 * Function which generically implements an element wise two-argument
	 * function application. Used to implement mathematical operations.
	 *
	 * @tparam Func is the type of the function that should be applied.
	 * @tparam Other is the type of the second vector.
	 * @param v1 is the vector containing the left-hand arguments.
	 * @param v2 is the vector containing the right-hand arguments.
	 * @param f is the function that should be applied to v1 and v2.
	 * @return a new vector instance.
	 */
	template <typename Func, typename Other>
	friend Inst map(const Self &v1, const Other &v2, Func f)
	{
		static_assert(Self::size() == Other::size(), "Vector size missmatch");

		Inst res;
		for (size_t i = 0; i < v1.size(); i++) {
			res[i] = f(v1[i], v2[i]);
		}
		return res;
	}

	template <typename Func>
	friend Inst map(const Self &v, Func f)
	{
		Inst res;
		for (size_t i = 0; i < v.size(); i++) {
			res[i] = f(v[i]);
		}
		return res;
	}

	template <typename Func, typename Other>
	friend void map_onto(Self &v1, const Other &v2, Func f)
	{
		static_assert(Self::size() == Other::size(), "Vector size missmatch");

		for (size_t i = 0; i < v1.size(); i++) {
			f(v1[i], v2[i]);
		}
	}

	template <typename Func>
	friend void map_onto(Self &v, Func f)
	{
		for (size_t i = 0; i < v.size(); i++) {
			f(v[i]);
		}
	}

	template <typename Func1, typename Func2, typename Other>
	friend bool fold_bool(const Self &v1, const Other &v2, bool init, Func1 f1,
	                      Func2 f2)
	{
		static_assert(Self::size() == Other::size(), "Vector size missmatch");

		bool res = init;
		for (size_t i = 0; i < v1.size(); i++) {
			res = f2(res, f1(v1[i], v2[i]));
		}
		return res;
	}

public:
	static constexpr size_t size() { return Size; }

	template <typename U = Impl, typename = typename std::enable_if<
	                                 !std::is_const<U>::value>::type>
	pointer begin()
	{
		return &mem()[0];
	}

	template <typename U = Impl, typename = typename std::enable_if<
	                                 !std::is_const<U>::value>::type>
	pointer end()
	{
		return &mem()[Size];
	}

	const_pointer begin() const { return &mem()[0]; }

	const_pointer end() const { return &mem()[Size]; }

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	void assign(const Other &v)
	{
		static_assert(Self::size() == Other::size(), "Vector size missmatch");
		std::copy(v.begin(), v.end(), begin());
	}

	void assign(value_type s)
	{
		map_onto(*this, [s](value_type &a) { return a = s; });
	}

	/**
	 * Returns a reference at the i-th element in the vector.
	 */
	value_type &operator[](size_t idx) { return mem()[idx]; }

	/**
	 * Returns a copy of the i-th element in the vector.
	 */
	value_type operator[](size_t idx) const { return mem()[idx]; }

	value_type sqrL2Norm() const
	{
		value_type res = 0;
		for (size_t i = 0; i < Size; i++) {
			res += (*this)[i] * (*this)[i];
		}
		return Size == 0 ? res : res * value_type(1.0 / double(Size));
	}

	value_type L2Norm() const { return std::sqrt(sqrL2Norm()); }

	friend std::ostream &operator<<(std::ostream &os, const Self &m)
	{
		for (size_t i = 0; i < Size; i++) {
			os << (i == 0 ? "" : ", ") << m[i];
		}
		return os;
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend void operator+=(Self &v1, const Other &v2)
	{
		map_onto(v1, v2, [](value_type &a, value_type b) { return a += b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend void operator-=(Self &v1, const Other &v2)
	{
		map_onto(v1, v2, [](value_type &a, value_type b) { return a -= b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend void operator*=(Self &v1, const Other &v2)
	{
		map_onto(v1, v2, [](value_type &a, value_type b) { return a *= b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend void operator/=(Self &v1, const Other &v2)
	{
		map_onto(v1, v2, [](value_type &a, value_type b) { return a /= b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend Inst operator+(const Self &v1, const Other &v2)
	{
		return map(v1, v2, [](value_type a, value_type b) { return a + b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend Inst operator-(const Self &v1, const Other &v2)
	{
		return map(v1, v2, [](value_type a, value_type b) { return a - b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend Inst operator*(const Self &v1, const Other &v2)
	{
		return map(v1, v2, [](value_type a, value_type b) { return a * b; });
	}

	template <typename Other, typename = typename std::enable_if<
	                              !std::is_scalar<Other>::value>::type>
	friend Inst operator/(const Self &v1, const Other &v2)
	{
		return map(v1, v2, [](value_type a, value_type b) { return a / b; });
	}

	friend void operator+=(Self &v, value_type s)
	{
		map_onto(v, [s](value_type &a) { return a += s; });
	}

	friend void operator-=(Self &v, value_type s)
	{
		map_onto(v, [s](value_type &a) { return a -= s; });
	}

	friend void operator*=(Self &v, value_type s)
	{
		map_onto(v, [s](value_type &a) { return a *= s; });
	}

	friend void operator/=(Self &v, value_type s)
	{
		map_onto(v, [s](value_type &a) { return a /= s; });
	}

	friend Inst operator+(const Self &v, value_type s)
	{
		return map(v, [s](value_type a) { return a + s; });
	}

	friend Inst operator-(const Self &v, value_type s)
	{
		return map(v, [s](value_type a) { return a - s; });
	}

	friend Inst operator*(value_type s, const Self &v)
	{
		return map(v, [s](value_type a) { return s * a; });
	}

	friend Inst operator*(const Self &v, value_type s)
	{
		return map(v, [s](value_type a) { return a * s; });
	}

	friend Inst operator/(const Self &v, value_type s)
	{
		return map(v, [s](value_type a) { return a / s; });
	}

	template <typename Other>
	friend bool operator==(const Self &v1, const Other &v2)
	{
		return fold_bool(v1, v2, true,
		                 [](value_type a, value_type b) { return a == b; },
		                 [](bool a, bool b) { return a && b; });
	}

	template <typename Other>
	friend bool operator!=(const Self &v1, const Other &v2)
	{
		return fold_bool(v1, v2, false,
		                 [](value_type a, value_type b) { return a != b; },
		                 [](bool a, bool b) { return a || b; });
	}

	/**
	 * Returns a view onto this vector of ViewSize elements, starting at the
	 * element indicated by ViewOffs.
	 */
	template <size_t ViewSize, size_t ViewOffs = 0>
	auto view()
	{
		return VectorView<T, ViewSize>(mem() + ViewOffs);
	}

	template <size_t ViewSize, size_t ViewOffs = 0>
	auto view() const
	{
		return VectorView<const T, ViewSize>(mem() + ViewOffs);
	}
};

/**
 * Vector class without own storage which implements a view at another vector.
 * The lifetime of a VectorView instance is coupled to the underlying vector.
 */
template <typename T, size_t Size>
class VectorView
    : public VectorMixin<VectorView<T, Size>, Vector<T, Size>, T, Size> {
private:
	friend class VectorMixin<VectorView<T, Size>, Vector<T, Size>, T, Size>;

	/**
	 * Pointer at the first storage element.
	 */
	T *m_mem;

	/**
	 * Returns the storage pointer.
	 */
	T *mem() { return m_mem; }

	/**
	 * Returns a const version of the storage pointer.
	 */
	const T *mem() const { return m_mem; }

public:
	/**
	 * Create a new VectorView instance with the first element pointing at the
	 * given memory-location.
	 *
	 * @param mem is the memory-location of the first element in the vector
	 * view.
	 */
	VectorView(T *mem) : m_mem(mem) {}
};

/**
 * Base class for custom fixed-size vector types. In contrast to VectorMixin and
 * VectorView this class provides its own storage space. Derive from this class
 * to provide a specialised vector class with its own member functions.
 */
template <typename Inst_, typename T, size_t Size>
class alignas(16) VectorBase
    : public VectorMixin<VectorBase<Inst_, T, Size>, Inst_, T, Size> {
private:
	friend class VectorMixin<VectorBase<Inst_, T, Size>, Inst_, T, Size>;

	/**
	 * Fixed-size array containing the vector elements.
	 */
	std::array<T, Size> m_arr;

	/**
	 * Provides access at the first element in the array.
	 */
	T *mem() { return &m_arr[0]; }

	/**
	 * Provides const access at the first element in the array.
	 */
	const T *mem() const { return &m_arr[0]; }

	/**
	 * Used internally to construct an array of names from the information
	 * returned from the info() method.
	 */
	template <size_t... Is>
	static constexpr std::array<const char *, Size> names(seq<Is...>)
	{
		return {{info<Is>().name...}};
	}

	/**
	 * Used internally to construct an array of units from the information
	 * returned from the info() method.
	 */
	template <size_t... Is>
	static constexpr std::array<const char *, Size> units(seq<Is...>)
	{
		return {{info<Is>().unit...}};
	}

	/**
	 * Used internally to construct an array of units from the information
	 * returned from the info() method.
	 */
	template <size_t... Is>
	static constexpr std::array<T, Size> scales(seq<Is...>)
	{
		return {{info<Is>().scale...}};
	}

protected:
	/**
	 * Small template type used for accessing element-specific information by
	 * the mechanism of function overloads.
	 */
	template <size_t I_>
	struct InfoAccessor {
		static constexpr size_t I = I_;
	};

public:
	using Self = VectorBase<Inst_, T, Size>;
	using Inst = Inst_;

	constexpr std::array<T, Size> as_array() const { return m_arr; }

	/**
	 * Default constructor. Default-initializes all vector elements.
	 */
	constexpr VectorBase() {}

	/**
	 * Constructor which allows to initialize the vector with an array
	 * containing the vector elements.
	 */
	constexpr VectorBase(const std::array<T, Size> &arr) : m_arr(arr) {}

	/**
	 * Internal method returning information about the I-th vector element. Do
	 * not use this function directly, use info instead.
	 *
	 * @tparam Accessor is a type which must provide an "I" constant specifying
	 * the index of the element for which the run-time or compile-time
	 * information should be returned.
	 */
	template <typename Accessor>
	static constexpr VectorElementInfo element_info(Accessor)
	{
		static_assert(Accessor::I < Size,
		              "Cannot access vector element info, index out of range");
		return VectorElementInfo();
	}

	/**
	 * Returns information about the I-th vector element.
	 *
	 * @tparam I is the index of the vector element for which information should
	 * be returned.
	 */
	template <size_t I>
	static constexpr VectorElementInfo info()
	{
		return Inst::element_info(InfoAccessor<I>());
	}

	/**
	 * Returns the names of the vector elements.
	 *
	 * @return an array of the size of the vector containing the name of each
	 * element.
	 */
	static constexpr std::array<const char *, Size> names()
	{
		return names(gen_seq<Size>());
	}

	/**
	 * Returns the units of the vector elements.
	 *
	 * @return an array of the size of the vector containing the physical unit
	 * of each element.
	 */
	static constexpr std::array<const char *, Size> units()
	{
		return units(gen_seq<Size>());
	}

	/**
	 * Returns the scales of the vector elements. The scale corresponds to the
	 * value each element has to be multiplied with to typically transform it
	 * into a value range where it can be represented without an exponent.
	 * Appart from visualisation purposes, the scales() method is used by the
	 * DormandPrince integrator to normalise the error values of each component.
	 *
	 * @return an instance of this vector where each element corresponds to the
	 * scalue the elements typically should be multiplied with.
	 */
	static constexpr Inst scales() { return Inst(scales(gen_seq<Size>())); }
};

/**
 * Simple specialization of VectorBase. Provides a mathematical vector with
 * fixed size. This class is also used to store intermediate results from vector
 * views.
 */
template <typename T, size_t Size>
class Vector final : public VectorBase<Vector<T, Size>, T, Size> {
public:
	using VectorBase<Vector, T, Size>::VectorBase;
};

/**
 * Class used to calculate the total dimensionality of a state vector.
 */
template <typename... Ts>
struct MultiVectorDim;

template <>
struct MultiVectorDim<> {
	static constexpr size_t calculate() { return 0; }
};

template <typename T0, typename... Ts>
struct MultiVectorDim<T0, Ts...> {
	static constexpr size_t calculate()
	{
		return T0::size() + MultiVectorDim<Ts...>::calculate();
	}
};

/**
 * Vector type combining multiple other vector types into one while still
 * allowing to call the introspection member functions such as names(), units(),
 * scales() and info<I>().
 */
template <typename Inst_, typename Vec0, typename... Vecs>
class MultiVector
    : public VectorBase<Inst_, typename Vec0::value_type,
                        MultiVectorDim<Vec0, Vecs...>::calculate()> {
public:
	using Base = VectorBase<Inst_, typename Vec0::value_type,
	                        MultiVectorDim<Vec0, Vecs...>::calculate()>;
	using Base::Base;

private:
	/**
	 * Helper structure used to make sure that the inner vector types do not
	 * posses any additional member variables, as this would interfere with the
	 * way the default constructors of the inner types are called.
	 */
	template <size_t Expected, size_t Actual>
	struct SizeCheck {
		/**
		 * Calculates the expected size taking the alignment into account.
		 */
		static constexpr size_t aligned_size(size_t s)
		{
			constexpr size_t ALIGN = 16;
			return s == 0 ? ALIGN : ((s + ALIGN - 1) / ALIGN) * ALIGN;
		}

		static_assert(
		    aligned_size(Expected) == aligned_size(Actual),
		    "Vectors used in conjunction with MultiVector may not possess any "
		    "additional data members!");
	};

	/**
	 * Recursion termination for the element_info introspection method. Simply
	 * calles the element_info method for the element with the index I - Offs
	 * in the vector Vec0.
	 */
	template <size_t I, size_t Offs>
	static constexpr VectorElementInfo element_info_impl()
	{
		return VectorElementInfo();
	}

	/**
	 * Actuall call to the info method of the internal vector type. Only enabled
	 * if the given index is valid.
	 */
	template <size_t I, typename V>
	static constexpr VectorElementInfo element_info_impl(
	    typename std::enable_if<(I < V::Size)>::type * = 0)
	{
		return V::template info<I>();
	}

	/**
	 * If the given index is not valid, return an empty VectorElementInfo
	 * object.
	 */
	template <size_t I, typename V>
	static constexpr VectorElementInfo element_info_impl(
	    typename std::enable_if<(I >= V::Size)>::type * = 0)
	{
		return VectorElementInfo();
	}

	/**
	 * Compile-time recursive method which obtains information about the I-th
	 * element by going over each vector the MultiVector class is constituted
	 * of.
	 */
	template <size_t I, size_t Offs, typename V0, typename... Vs>
	static constexpr VectorElementInfo element_info_impl()
	{
		return I < V0::Size + Offs
		           ? element_info_impl<I - Offs, V0>()
		           : element_info_impl<I, Offs + V0::Size, Vs...>();
	}

	/**
	 * Recursion termination of the compile-time recursive implementation of
	 * ctor_impl().
	 */
	template <size_t Offs>
	void ctor_impl(typename Base::pointer)
	{
	}

	/**
	 * Compile-time recursive implementation of the default constructor.
	 * Constructs the inner vectors at the corresponding memory addresses.
	 */
	template <size_t Offs, typename V0, typename... Vs>
	void ctor_impl(typename Base::pointer tar)
	{
		// Call the constructor at the target memory address.
		// Note: This is a very dirty hack which may cause strange bugs as soon
		// as child classes have additional data members. However, any other
		// means of calling the V0 constructor caused a factor two slowdown in
		// the unit tests which is not acceptable.
		SizeCheck<sizeof(V0), V0::Size * sizeof(typename Base::value_type)>();
		new (tar + Offs) V0;
		ctor_impl<Offs + V0::Size, Vs...>(tar);
	}

public:
	/**
	 * Default constructor: recursively call the inner default constructors.
	 */
	MultiVector() { ctor_impl<0, Vec0, Vecs...>(Base::begin()); }

	/**
	 * Internal method returning information about the I-th vector element. Do
	 * not use this function directly.
	 *
	 * @tparam Accessor is a type which must provide an "I" constant specifying
	 * the index of the element for which the run-time or compile-time
	 * information should be returned.
	 */
	template <typename Accessor>
	static constexpr VectorElementInfo element_info(Accessor)
	{
		static_assert(Accessor::I < Base::Size,
		              "Cannot access vector element info, index out of range");
		return element_info_impl<Accessor::I, 0, Vec0, Vecs...>();
	}
};

#define NAMED_VECTOR_ELEMENT(NAME, IDX)                                \
	static constexpr size_t idx_##NAME = IDX;                          \
	static constexpr VectorElementInfo element_info(InfoAccessor<IDX>) \
	{                                                                  \
		return VectorElementInfo(#NAME);                               \
	}                                                                  \
	Inst &NAME(value_type x)                                           \
	{                                                                  \
		(*this)[IDX] = x;                                              \
		return static_cast<Inst &>(*this);                             \
	}                                                                  \
	value_type &NAME() { return (*this)[IDX]; }                        \
	value_type NAME() const { return (*this)[IDX]; }

#define TYPED_VECTOR_ELEMENT(NAME, IDX, TYPE)                          \
	static constexpr size_t idx_##NAME = IDX;                          \
	static constexpr VectorElementInfo element_info(InfoAccessor<IDX>) \
	{                                                                  \
		return VectorElementInfo(#NAME, TYPE::unit(), TYPE::scale());  \
	}                                                                  \
	Inst &NAME(TYPE x)                                                 \
	{                                                                  \
		(*this)[IDX] = x.v();                                          \
		return static_cast<Inst &>(*this);                             \
	}                                                                  \
	TYPE NAME() const { return TYPE((*this)[IDX]); }
}

#endif /* CINDER_COMMON_VECTOR_HPP */
