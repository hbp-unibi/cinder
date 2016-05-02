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

#include "gtest/gtest.h"

#include <cinder/common/vector.hpp>
#include <cinder/common/types.hpp>

namespace cinder {

TEST(vector, view)
{
	Vector<int, 8> vec({0, 1, 2, 3, 4, 5, 6, 7});
	auto view = vec.view<4, 3>();

	for (int i = 0; i < int(view.size()); i++) {
		EXPECT_EQ(i + 3, view[i]);
	}

	auto view2 = view.view<2, 2>();
	for (int i = 0; i < int(view2.size()); i++) {
		EXPECT_EQ(i + 5, view2[i]);
	}
}

TEST(vector, iterator)
{
	Vector<int, 8> vec({0, 1, 2, 3, 4, 5, 6, 7});
	int i = 0;
	for (int j : vec) {
		EXPECT_EQ(i++, j);
	}
	EXPECT_EQ(8, i);

	EXPECT_EQ(8, vec.end() - vec.begin());
}

TEST(vector, writeable_iterator)
{
	Vector<int, 8> vec({0, 1, 2, 3, 4, 5, 6, 7});

	*(vec.begin()) = 5;
	EXPECT_EQ(5, vec[0]);
}

TEST(vector, writeable_view_iterator)
{
	Vector<int, 8> vec({0, 1, 2, 3, 4, 5, 6, 7});

	*(vec.view<4, 3>().begin()) = 5;
	EXPECT_EQ(5, vec[3]);
}

TEST(vector, scalar_arithmetic)
{
	Vector<int, 4> vec({0, 1, 2, 3});
	auto view = vec.view<2, 1>();

	{
		auto v = vec + 5;
		for (int i = 0; i < int(v.size()); i++) {
			EXPECT_EQ(i + 5, v[i]);
		}
	}

	{
		auto v = view + 5;
		for (int i = 0; i < int(v.size()); i++) {
			EXPECT_EQ(i + 6, v[i]);
		}
	}

	{
		auto v = vec - 5;
		for (int i = 0; i < int(v.size()); i++) {
			EXPECT_EQ(i - 5, v[i]);
		}
	}

	{
		auto v = view - 5;
		for (int i = 0; i < int(v.size()); i++) {
			EXPECT_EQ(i - 4, v[i]);
		}
	}

	{
		auto v = vec * 2;
		for (int i = 0; i < int(v.size()); i++) {
			EXPECT_EQ(2 * i, v[i]);
		}
	}

	{
		auto v = view * 2;
		for (int i = 0; i < int(v.size()); i++) {
			EXPECT_EQ(2 * (i + 1), v[i]);
		}
	}

	{
		auto v = vec / 2;
		for (int i = 0; i < int(v.size()); i++) {
			EXPECT_EQ(i / 2, v[i]);
		}
	}

	{
		auto v = view / 2;
		for (int i = 0; i < int(v.size()); i++) {
			EXPECT_EQ((i + 1) / 2, v[i]);
		}
	}
}

TEST(vector, vector_arithmetic)
{
	Vector<int, 4> vec1({0, 1, 2, 3});
	Vector<int, 4> vec2 = vec1 + 1;

	{
		auto v = vec1 + vec2;
		for (size_t i = 0; i < v.size(); i++) {
			EXPECT_EQ(vec1[i] + vec2[i], v[i]);
		}
	}

	{
		auto v = vec1 * vec2;
		for (size_t i = 0; i < v.size(); i++) {
			EXPECT_EQ(vec1[i] * vec2[i], v[i]);
		}
	}

	{
		auto v = vec1 / vec2;
		for (size_t i = 0; i < v.size(); i++) {
			EXPECT_EQ(vec1[i] / vec2[i], v[i]);
		}
	}

	{
		auto v = vec1 - vec2;
		for (size_t i = 0; i < v.size(); i++) {
			EXPECT_EQ(vec1[i] - vec2[i], v[i]);
		}
	}
}

TEST(vector, scalar_inplace_arithmetic)
{
	Vector<size_t, 4> vec({0, 1, 2, 3});
	auto view = vec.view<2, 1>();

	for (size_t i = 0; i < vec.size(); i++) {
		EXPECT_EQ(i, vec[i]);
	}
	for (size_t i = 0; i < view.size(); i++) {
		EXPECT_EQ(i + 1, view[i]);
	}

	vec += 5;
	for (size_t i = 0; i < vec.size(); i++) {
		EXPECT_EQ(i + 5, vec[i]);
	}
	for (size_t i = 0; i < view.size(); i++) {
		EXPECT_EQ(i + 6, view[i]);
	}

	vec *= 2;
	for (size_t i = 0; i < vec.size(); i++) {
		EXPECT_EQ(2 * (i + 5), vec[i]);
	}
	for (size_t i = 0; i < view.size(); i++) {
		EXPECT_EQ(2 * (i + 6), view[i]);
	}

	vec /= 2;
	for (size_t i = 0; i < vec.size(); i++) {
		EXPECT_EQ(i + 5, vec[i]);
	}
	for (size_t i = 0; i < view.size(); i++) {
		EXPECT_EQ(i + 6, view[i]);
	}

	vec -= 2;
	for (size_t i = 0; i < vec.size(); i++) {
		EXPECT_EQ(i + 3, vec[i]);
	}
	for (size_t i = 0; i < view.size(); i++) {
		EXPECT_EQ(i + 4, view[i]);
	}

	// Test operations on the view
	view += 2;
	for (size_t i = 0; i < vec.size(); i++) {
		EXPECT_EQ((i < 1 || i > 2) ? i + 3 : i + 5, vec[i]);
	}
	for (size_t i = 0; i < view.size(); i++) {
		EXPECT_EQ(i + 6, view[i]);
	}

	view -= 2;
	for (size_t i = 0; i < vec.size(); i++) {
		EXPECT_EQ(i + 3, vec[i]);
	}
	for (size_t i = 0; i < view.size(); i++) {
		EXPECT_EQ(i + 4, view[i]);
	}

	view *= 2;
	for (size_t i = 0; i < vec.size(); i++) {
		EXPECT_EQ((i < 1 || i > 2) ? i + 3 : 2 * (i + 3), vec[i]);
	}
	for (size_t i = 0; i < view.size(); i++) {
		EXPECT_EQ(2 * (i + 4), view[i]);
	}

	view /= 2;
	for (size_t i = 0; i < vec.size(); i++) {
		EXPECT_EQ(i + 3, vec[i]);
	}
	for (size_t i = 0; i < view.size(); i++) {
		EXPECT_EQ(i + 4, view[i]);
	}
}

TEST(vector, vector_inplace_arithmetic)
{
	Vector<int, 4> vec1({0, 1, 2, 3});
	Vector<int, 4> vec2 = vec1 + 1;

	{
		auto v = vec1;
		v += vec2;
		for (size_t i = 0; i < v.size(); i++) {
			EXPECT_EQ(vec1[i] + vec2[i], v[i]);
		}
	}

	{
		auto v = vec1;
		v *= vec2;
		for (size_t i = 0; i < v.size(); i++) {
			EXPECT_EQ(vec1[i] * vec2[i], v[i]);
		}
	}

	{
		auto v = vec1;
		v /= vec2;
		for (size_t i = 0; i < v.size(); i++) {
			EXPECT_EQ(vec1[i] / vec2[i], v[i]);
		}
	}

	{
		auto v = vec1;
		v -= vec2;
		for (size_t i = 0; i < v.size(); i++) {
			EXPECT_EQ(vec1[i] - vec2[i], v[i]);
		}
	}
}

TEST(vector, equality)
{
	Vector<int, 3> v1({1, 2, 3});
	Vector<int, 3> v2({1, 2, 3});
	Vector<int, 3> v3({2, 2, 3});
	Vector<int, 3> v4({1, 3, 3});
	Vector<int, 3> v5({1, 2, 4});

	EXPECT_TRUE(v1 == v2);
	EXPECT_FALSE(v1 != v2);

	EXPECT_TRUE(v1 != v3);
	EXPECT_FALSE(v1 == v3);

	EXPECT_TRUE(v1 != v4);
	EXPECT_FALSE(v1 == v4);

	EXPECT_TRUE(v1 != v5);
	EXPECT_FALSE(v1 == v5);
}

namespace {
struct TestVector : public VectorBase<TestVector, Real, 6> {
	using VectorBase<TestVector, Real, 6>::VectorBase;

	TestVector() { a(42.0); }

	NAMED_VECTOR_ELEMENT(a, 0);
	TYPED_VECTOR_ELEMENT(u, 1, Voltage);
	TYPED_VECTOR_ELEMENT(i, 2, Current);
	TYPED_VECTOR_ELEMENT(c, 3, Capacitance);
	TYPED_VECTOR_ELEMENT(g, 4, Conductance);
	TYPED_VECTOR_ELEMENT(t, 5, RealTime);
};

struct TestVector2 : public VectorBase<TestVector2, Real, 3> {
	using VectorBase<TestVector2, Real, 3>::VectorBase;

	TestVector2() { b(84.0); }

	NAMED_VECTOR_ELEMENT(b, 0);
	TYPED_VECTOR_ELEMENT(c, 1, Voltage);
	TYPED_VECTOR_ELEMENT(d, 2, Current);
};

struct TestMultiVector
    : public MultiVector<TestMultiVector, TestVector, TestVector2> {
	using MultiVector<TestMultiVector, TestVector, TestVector2>::MultiVector;
};

static_assert(std::is_same<TestMultiVector::value_type, Real>::value,
              "Invalid multi vector value type");
static_assert(TestMultiVector::Size == 9, "Invalid multi vector size");
}

TEST(vector, named_elements)
{
	EXPECT_EQ("a", TestVector::info<0>().name);
	EXPECT_EQ("u", TestVector::info<1>().name);
	EXPECT_EQ("i", TestVector::info<2>().name);
	EXPECT_EQ("c", TestVector::info<3>().name);
	EXPECT_EQ("g", TestVector::info<4>().name);
	EXPECT_EQ("t", TestVector::info<5>().name);

	EXPECT_EQ("", TestVector::info<0>().unit);
	EXPECT_EQ("V", TestVector::info<1>().unit);
	EXPECT_EQ("A", TestVector::info<2>().unit);
	EXPECT_EQ("F", TestVector::info<3>().unit);
	EXPECT_EQ("S", TestVector::info<4>().unit);
	EXPECT_EQ("s", TestVector::info<5>().unit);

	EXPECT_EQ(1e0, TestVector::info<0>().scale);
	EXPECT_EQ(1e3, TestVector::info<1>().scale);
	EXPECT_EQ(1e9, TestVector::info<2>().scale);
	EXPECT_EQ(1e9, TestVector::info<3>().scale);
	EXPECT_EQ(1e6, TestVector::info<4>().scale);
	EXPECT_EQ(1e3, TestVector::info<5>().scale);

	std::array<const char *, 6> expected_names = {"a", "u", "i", "c", "g", "t"};
	std::array<const char *, 6> expected_units = {"", "V", "A", "F", "S", "s"};
	TestVector expected_scales = {{1e0, 1e3, 1e9, 1e9, 1e6, 1e3}};
	std::array<VectorElementInfo, 6> expected_infos = {{
	    {"a", "", 1e0},  {"u", "V", 1e3}, {"i", "A", 1e9},
	    {"c", "F", 1e9}, {"g", "S", 1e6}, {"t", "s", 1e3}}};

	EXPECT_EQ(expected_names, TestVector::names());
	EXPECT_EQ(expected_units, TestVector::units());
	EXPECT_EQ(expected_scales, TestVector::scales());
	EXPECT_EQ(expected_infos, TestVector::infos());
}

TEST(multi_vector, named_elements)
{
	EXPECT_EQ("a", TestMultiVector::info<0>().name);
	EXPECT_EQ("u", TestMultiVector::info<1>().name);
	EXPECT_EQ("i", TestMultiVector::info<2>().name);
	EXPECT_EQ("c", TestMultiVector::info<3>().name);
	EXPECT_EQ("g", TestMultiVector::info<4>().name);
	EXPECT_EQ("t", TestMultiVector::info<5>().name);
	EXPECT_EQ("b", TestMultiVector::info<6>().name);
	EXPECT_EQ("c", TestMultiVector::info<7>().name);
	EXPECT_EQ("d", TestMultiVector::info<8>().name);

	EXPECT_EQ("", TestMultiVector::info<0>().unit);
	EXPECT_EQ("V", TestMultiVector::info<1>().unit);
	EXPECT_EQ("A", TestMultiVector::info<2>().unit);
	EXPECT_EQ("F", TestMultiVector::info<3>().unit);
	EXPECT_EQ("S", TestMultiVector::info<4>().unit);
	EXPECT_EQ("s", TestMultiVector::info<5>().unit);
	EXPECT_EQ("", TestMultiVector::info<6>().unit);
	EXPECT_EQ("V", TestMultiVector::info<7>().unit);
	EXPECT_EQ("A", TestMultiVector::info<8>().unit);

	EXPECT_EQ(1e0, TestMultiVector::info<0>().scale);
	EXPECT_EQ(1e3, TestMultiVector::info<1>().scale);
	EXPECT_EQ(1e9, TestMultiVector::info<2>().scale);
	EXPECT_EQ(1e9, TestMultiVector::info<3>().scale);
	EXPECT_EQ(1e6, TestMultiVector::info<4>().scale);
	EXPECT_EQ(1e3, TestMultiVector::info<5>().scale);
	EXPECT_EQ(1e0, TestMultiVector::info<6>().scale);
	EXPECT_EQ(1e3, TestMultiVector::info<7>().scale);
	EXPECT_EQ(1e9, TestMultiVector::info<8>().scale);

	std::array<const char *, 9> expected_names = {"a", "u", "i", "c", "g",
	                                              "t", "b", "c", "d"};
	std::array<const char *, 9> expected_units = {"",  "V", "A", "F", "S",
	                                              "s", "",  "V", "A"};
	TestMultiVector expected_scales = {
	    {1e0, 1e3, 1e9, 1e9, 1e6, 1e3, 1e0, 1e3, 1e9}};
	std::array<VectorElementInfo, 9> expected_infos = {{
	    {"a", "", 1e0},  {"u", "V", 1e3}, {"i", "A", 1e9},
	    {"c", "F", 1e9}, {"g", "S", 1e6}, {"t", "s", 1e3},
	    {"b", "", 1e0},  {"c", "V", 1e3}, {"d", "A", 1e9}}};

	EXPECT_EQ(expected_names, TestMultiVector::names());
	EXPECT_EQ(expected_units, TestMultiVector::units());
	EXPECT_EQ(expected_scales, TestMultiVector::scales());
	EXPECT_EQ(expected_infos, TestMultiVector::infos());

	TestMultiVector mv;

	EXPECT_EQ(42.0_R, mv[0]);
	EXPECT_EQ(84.0_R, mv[6]);
}

TEST(multi_vector, get)
{
	TestMultiVector mv({{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}}, {{7.0, 8.0, 9.0}});

	TestVector v1 = mv.get<0>();
	TestVector2 v2 = mv.get<1>();

	EXPECT_EQ(1.0_R, v1[0]);
	EXPECT_EQ(2.0_R, v1[1]);
	EXPECT_EQ(3.0_R, v1[2]);
	EXPECT_EQ(4.0_R, v1[3]);
	EXPECT_EQ(5.0_R, v1[4]);
	EXPECT_EQ(6.0_R, v1[5]);

	EXPECT_EQ(7.0_R, v2[0]);
	EXPECT_EQ(8.0_R, v2[1]);
	EXPECT_EQ(9.0_R, v2[2]);
}
}

