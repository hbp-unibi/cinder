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
	for (int j: vec) {
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
}

