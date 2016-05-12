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
 * @file spikes.hpp
 *
 * Contains methods which allow to construct and concatenate spike trains.
 *
 * @author Andreas Stöckel
 */

#include <algorithm>
#include <vector>
#include <random>

#include <cinder/common/types.hpp>

namespace cinder {
/**
 * Inserts spikes with a constant interval into the given spike train.
 *
 * @param spike_train is the spike train into which the spikes should be
 * inserted.
 * @param t_start is the start time.
 * @param t_end is the end time.
 * @param interval is the interval between two spikes.
 * @return a vector containing the spike times.
 */
static inline std::vector<Time> &constant_interval_spike_train(
    std::vector<Time> &spike_train, Time t_start, Time t_end, Time interval)
{
	for (Time t = t_start; t < t_end; t += interval) {
		spike_train.emplace_back(t);
	}
	return spike_train;
}

/**
 * Creates a spike train with spikes starting at the given time, ending at the
 * given time and with the specified inter-spike-interval.
 *
 * @param t_start is the start time.
 * @param t_end is the end time.
 * @param interval is the interval between two spikes.
 * @return a vector containing the spike times.
 */
static inline std::vector<Time> constant_interval_spike_train(Time t_start,
                                                              Time t_end,
                                                              Time interval)
{
	std::vector<Time> spike_train;
	constant_interval_spike_train(spike_train, t_start, t_end, interval);
	return spike_train;
}

/**
 * Inserts n spikes with a constant interval into the given spike train.
 *
 * @param spike_train is the spike train into which the spikes should be
 * inserted.
 * @param t_start is the start time.
 * @param n is the number of spikes that should be inserted.
 * @param interval is the interval between two spikes.
 * @return a vector containing the spike times.
 */
static inline std::vector<Time> &constant_interval_spike_train(
    std::vector<Time> &spike_train, Time t_offs, size_t n, Time interval)
{
	Time t = t_offs;
	for (size_t i = 0; i < n; i++, t += interval) {
		spike_train.emplace_back(t);
	}
	return spike_train;
}

/**
 * Creates a spike train with n spikes starting at the given time with the
 * specified inter-spike-interval.
 *
 * @param t_start is the start time.
 * @param t_end is the end time.
 * @param interval is the interval between two spikes.
 * @return a vector containing the spike times.
 */
static inline std::vector<Time> constant_interval_spike_train(Time t_offs,
                                                              size_t n,
                                                              Time interval)
{
	std::vector<Time> spike_train;
	constant_interval_spike_train(spike_train, t_offs, n, interval);
	return spike_train;
}

/**
 * Creates a spike train with spikes occuring with the given frequency.
 *
 * @param t_start is the start time.
 * @param t_end is the end time.
 * @param frequency is the spike frequency in hertz (per second).
 * @return a vector containing the spike times.
 */
static inline std::vector<Time> &constant_frequency_spike_train(
    std::vector<Time> &spike_train, Time t_start, Time t_end, float frequency)
{
	return constant_interval_spike_train(spike_train, t_start, t_end,
	                                     Time::sec(1.0 / frequency));
}

/**
 * Creates a spike train with spikes occuring with the given frequency.
 *
 * @param t_start is the start time.
 * @param t_end is the end time.
 * @param frequency is the spike frequency in hertz (per second).
 * @return a vector containing the spike times.
 */
static inline std::vector<Time> constant_frequency_spike_train(Time t_start,
                                                               Time t_end,
                                                               float frequency)
{
	return constant_interval_spike_train(t_start, t_end,
	                                     Time::sec(1.0 / frequency));
}

/**
 * Creates a spike train with n spikes occuring with the given frequency.
 *
 * @param t_offs is the start time.
 * @param n is the number of spikes that should be produced.
 * @param frequency is the spike frequency in hertz (per second).
 * @return a vector containing the spike times.
 */
static inline std::vector<Time> &constant_frequency_spike_train(
    std::vector<Time> &spike_train, Time t_offs, size_t n, float frequency)
{
	return constant_interval_spike_train(spike_train, t_offs, n,
	                                     Time::sec(1.0 / frequency));
}

/**
 * Creates a spike train with n spikes occuring with the given frequency.
 *
 * @param t_offs is the start time.
 * @param n is the number of spikes that should be produced.
 * @param frequency is the spike frequency in hertz (per second).
 * @return a vector containing the spike times.
 */
static inline std::vector<Time> constant_frequency_spike_train(Time t_offs,
                                                               size_t n,
                                                               float frequency)
{
	return constant_interval_spike_train(t_offs, n, Time::sec(1.0 / frequency));
}

/**
 * Adds gaussian jitter to the given spike train.
 *
 * @param spike_train is the spike train to which the jitter should be added.
 * @param sigma is the standard deviation of the Gaussian jitter that should be
 * added to the spike train.
 * @param seed if -1 a random seed is chosen, otherwise the given seed is used.
 * @return a reference to the given spike train.
 */
static inline std::vector<Time> &add_gaussian_jitter_to_spike_train(
    std::vector<Time> &spike_train, Time sigma, int seed = -1)
{
	std::default_random_engine gen(seed == -1 ? std::random_device()() : seed);
	std::normal_distribution<double> dist(0.0, sigma.sec());
	for (size_t i = 0; i < spike_train.size(); i++) {
		spike_train[i] += Time::sec(dist(gen));
	}
	return spike_train;
}

/**
 * Normalises the given spike train without creating a copy. Shifts the first
 * spike time to the given offset and sorts the spikes.
 *
 * @param spike_train is the spike_train that should be normalised.
 * @return a reference to the spike train passed as first parameter.
 */
static inline std::vector<Time> &normalise_spike_train(
    std::vector<Time> &spike_train, Time t_offs = 0_s)
{
	std::sort(spike_train.begin(), spike_train.end());
	for (size_t i = 0; i < spike_train.size(); i++) {
		spike_train[i] = spike_train[i] - spike_train[0] + t_offs;
	}
	return spike_train;
}
}
