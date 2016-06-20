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
 * @file noise_source.hpp
 *
 * Implementation of current sources which inject noise into the neuron model.
 *
 * @author Andreas Stöckel
 */

#ifndef CINDER_MODELS_NOISE_SOURCE_HPP
#define CINDER_MODELS_NOISE_SOURCE_HPP

#include <array>
#include <cmath>
#include <random>

#include <cinder/models/current_source.hpp>

namespace cinder {

/**
 * GaussianNoiseSourceState contains the state of the gaussian current source,
 * which consists of the current that was drawn at the last sampling time step.
 */
struct GaussianNoiseSourceState
    : public VectorBase<GaussianNoiseSourceState, Real, 1> {
	using VectorBase<GaussianNoiseSourceState, Real, 1>::VectorBase;

	TYPED_VECTOR_ELEMENT(i_noise, 0, Current);
};

/**
 * GaussianNoiseSourceParameters contains the parameters describing the behavior
 * of the GaussianNoiseSource. There are three parameters: i_stddev() which
 * defines the amplitude of the noise and seed() which is used as initial value
 * of the internally used random number generator. Furthermore, t_delta()
 * specifies the time interval at which the Gaussian noise is sampled, which
 * implicitly limits the bandwidth of the noise signal.
 */
struct GaussianNoiseSourceParameters
    : public VectorBase<GaussianNoiseSourceParameters, Real, 3> {
	using VectorBase<GaussianNoiseSourceParameters, Real, 3>::VectorBase;

	/**
	 * Constructor setting the default parameters.
	 */
	GaussianNoiseSourceParameters()
	{
		i_stddev(1_pA);
		t_delta(1_ms);
		seed(156484.0);
	}

	TYPED_VECTOR_ELEMENT(i_stddev, 0, Current);
	TYPED_VECTOR_ELEMENT(t_delta, 1, RealTime);
	NAMED_VECTOR_ELEMENT(seed, 2);
};

/**
 * A current source producing Gaussian white noise with the given standard
 * deviation. Note that this current source does not implement an autonomous
 * differential equation -- the "current" method implicitly depends on the
 * current time.
 */
class GaussianNoiseSource
    : public CurrentSourceBase<GaussianNoiseSourceState,
                               GaussianNoiseSourceParameters> {
private:
	std::default_random_engine m_gen;
	Time m_next_t;

public:
	using Base = CurrentSourceBase<GaussianNoiseSourceState,
	                               GaussianNoiseSourceParameters>;
	using Base::p;
	using Base::Base;

	/**
	 * Initializes the internal random number generator.
	 */
	template <typename State, typename System>
	void init(Time t, const State &, const System &)
	{
		m_gen.seed(p().seed());
		m_next_t = t;
	}

	/**
	 * Returns the position of the next discontinuity in the differential
	 * equation.
	 */
	Time next_discontinuity(Time) const { return m_next_t; }

	/**
	 * Called whenever the discontinuity is reached. Samples a new random value.
	 */
	template <typename State2, typename System>
	void handle_discontinuity(Time t, State2 &s, System &)
	{
		if (t >= m_next_t) {
			const RealTime dt = p().t_delta();
			const Real stddev = p().i_stddev() * std::sqrt(dt);
			std::normal_distribution<Real> dist(0.0, stddev);
			s[0] = dist(m_gen);
			m_next_t += Time::sec(dt);
		}
	}
};
}

#endif /* CINDER_MODELS_NOISE_SOURCE_HPP */
