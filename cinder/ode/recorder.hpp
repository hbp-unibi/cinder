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

#ifndef CINDER_ODE_RECORDER_HPP
#define CINDER_ODE_RECORDER_HPP

#include <iostream>

#include <cinder/common/time.hpp>

namespace cinder {
struct NullRecorder {
public:
	template <typename State, typename System>
	void record(Time, const State &, const System &)
	{
		// Do nothing here
	}
};

struct CSVRecorder {
private:
	std::ostream &m_os;
	Time m_last_time;
	Time m_min_delta;

public:
	CSVRecorder(std::ostream &os, Time min_delta = 0.1_ms): m_os(os), m_last_time(MIN_TIME), m_min_delta(min_delta) {}

	template <typename State, typename System>
	void record(Time t, const State &s, const System &)
	{
		// In case record is called multiple times in a row with the same
		// timestamp, just increment t by one granule and continue. Otherwise,
		// if less than m_min_delta has passed since the last recording, abort.
		if (m_last_time >= t) {
			t = m_last_time + Time(1);
		} else if (m_last_time + m_min_delta > t) {
			return;
		}
		m_os << t << ", " << s << '\n';
		m_last_time = t;
	}
};

}

#endif /* CINDER_MODELS_NEURONS_LIF_HPP */
