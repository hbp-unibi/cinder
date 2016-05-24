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

/**
 * Generates n input bursts with the given inter-spike-interval and the given
 * offset and spike jitter.
 *
 * @param spike_train is the spike train into which the burst should be
 * inserted.
 * @param burst_count is the number of bursts that should be generated.
 * @param burst_size is the number of spikes in one burst.
 * @param isi is the inter-spike interval within one burst
 * @param sigma_spike is the variance of the single spike times.
 * @param sigma_offs is the variance of the spike onset.
 * @param t_offs is the offset.
 * @param random if true, entirely random bursts are generated, otherwise a
 * representative burst is created.
 * @param seed if -1 a random seed is chosen, otherwise the given seed is used.
 */
static inline std::vector<Time> &generate_bursts(
    std::vector<Time> &spike_train, size_t burst_count, size_t burst_size,
    Time isi, Time sigma_spike, Time sigma_offs = 0_s, Time t_offs = 0_s,
    bool random = true, int seed = -1)
{
	// Calculate the deltaT for equidistant distribution
	const Time delta_t_eqn = Time(2 * (sigma_spike.t + sigma_offs.t) /
	                              std::max<size_t>(1, burst_count));

	// Random number generator
	std::default_random_engine gen(seed == -1 ? std::random_device()() : seed);
	std::normal_distribution<double> dist_spike(0.0, sigma_spike.sec());
	std::normal_distribution<double> dist_offs(0.0, sigma_offs.sec());

	// Iterate over all bursts and assemble the spikes
	for (size_t i = 0; i < burst_count; i++) {
		const Time offs = t_offs + (random ? Time::sec(dist_offs(gen)) : 0_s);
		for (size_t j = 0; j < burst_size; j++) {
			const Time t = Time(offs.t + isi.t * j);
			spike_train.emplace_back(t + (random
			                                  ? Time::sec(dist_spike(gen))
			                                  : Time(delta_t_eqn.t * i)));
		}
	}
	return spike_train;
}

/**
 * Generates n input bursts with the given inter-spike-interval and the given
 * offset and spike jitter and returns them.
 *
 * @param spike_train is the spike train into which the burst should be
 * inserted.
 * @param burst_count is the number of bursts that should be generated.
 * @param burst_size is the number of spikes in one burst.
 * @param isi is the inter-spike interval within one burst
 * @param sigma_spike is the variance of the single spike times.
 * @param sigma_offs is the variance of the spike onset.
 * @param t_offs is the offset.
 * @param random if true, entirely random bursts are generated, otherwise a
 * representative burst is created.
 * @param seed if -1 a random seed is chosen, otherwise the given seed is used.
 * @return the generated spike train
 */
static inline std::vector<Time> generate_bursts(
    size_t burst_count, size_t burst_size, Time isi, Time sigma_spike,
    Time sigma_offs = 0_s, Time t_offs = 0_s, bool random = true, int seed = -1)
{
	std::vector<Time> spike_train;
	generate_bursts(spike_train, burst_count, burst_size, isi, sigma_spike,
	                sigma_offs, t_offs, random, seed);
	return spike_train;
}
}
