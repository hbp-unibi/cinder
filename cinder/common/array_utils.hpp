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
 * @file lif.hpp
 *
 * Implementation of the simple linear integrate and fire neuron model.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_COMMON_ARRAY_UTILS_HPP
#define CINDER_COMMON_ARRAY_UTILS_HPP

#include <array>
#include <cstddef>

namespace cinder {

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

/**
 * Concatenates two arrays at compile time.
 */
template <typename T, size_t N1, size_t... I1, size_t N2, size_t... I2>
constexpr std::array<T, N1 + N2> concat(const std::array<T, N1> &a1,
                                        const std::array<T, N2> &a2, seq<I1...>,
                                        seq<I2...>)
{
	return {a1[I1]..., a2[I2]...};
}

/**
 * Concatenates two arrays at compile time.
 */
template <typename T, size_t N1, size_t N2>
constexpr std::array<T, N1 + N2> concat(const std::array<T, N1> &a1,
                                        const std::array<T, N2> &a2)
{
	return concat(a1, a2, gen_seq<N1>{}, gen_seq<N2>{});
}
}

#endif /* CINDER_COMMON_ARRAY_UTILS_HPP */
