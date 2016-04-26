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
 * @file recorder.hpp
 *
 * Contains some basic recorder classes. Recorders can be used to record the
 * state of the differential equation being integrated over time. To this end,
 * recorders provide a record() method which is called after each integration
 * timestep and after a discontinuity has been processed.
 *
 * Multiple Recorder instances can be combined using the MultiRecorder class and
 * the corresponding make_multi_recorder method.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_ODE_RECORDER_HPP
#define CINDER_ODE_RECORDER_HPP

#include <iostream>
#include <limits>

#include <cinder/common/time.hpp>
#include <cinder/common/types.hpp>

namespace cinder {
/**
 * The NullRecorder is the most basic recorder class which just does nothing and
 * exhibits no overhead at all.
 */
struct NullRecorder {
	/**
	 * The record method is called whenever a time slice was processed by the
	 * integrator. As its name suggests, the NullRecorder will just ignore the
	 * call to this method an do nothing.
	 */
	template <typename State, typename System>
	static void record(Time, const State &, const System &)
	{
		// Do nothing here
	}
};

/**
 * The CSVRecorder class can be used to write all neuron state variables over
 * time as CSV to a stream.
 */
class CSVRecorder {
private:
	/**
	 * Reference at the output stream to which the data should be written.
	 */
	std::ostream &m_os;

	/**
	 * Last time a datum was written to the output stream.
	 */
	Time m_last_time;

	/**
	 * User defined minimum delay between writing to points to the output
	 * stream.
	 */
	Time m_min_delta;

public:
	/**
	 * Constructor of the CSVRecorder class.
	 *
	 * @param os is the output stream to which the CSV should be written.
	 * @param min_delta is the minimum time difference between two samples in
	 * the output file. Note that min_delta is ignored in case samples are
	 * recorded with the same timestamp. This happens if a discontinuity is
	 * processed which should be recorded.
	 */
	CSVRecorder(std::ostream &os, Time min_delta = 0.1_ms)
	    : m_os(os), m_last_time(MIN_TIME), m_min_delta(min_delta)
	{
	}

	/**
	 * Called whenever the integrator processed a time slice. The CSVRecorder
	 * class will check whether at least m_min_delta (given in the constructor)
	 * has passed since something was written to the output stream. If yes, the
	 * current time stamp and the entire ODE state vector are written to the
	 * output.
	 */
	template <typename State, typename System>
	void record(Time t, const State &s, const System &)
	{
		// In case record is called multiple times in a row with the same
		// timestamp, just increment t by one granule and continue. Otherwise,
		// if less than m_min_delta has passed since the last recording, abort.
		if (m_last_time >= t) {
			t = m_last_time + Time(1);
		}
		else if (m_last_time + m_min_delta > t) {
			return;
		}
		m_os << t << ", " << s << '\n';
		m_last_time = t;
	}
};

/**
 * The MaximumMembranePotentialRecorder class records the maximum voltage
 * reached in a neuron simulation.
 */
struct MaximumMembranePotentialRecorder {
	/**
	 * Time at which the maximum occurred.
	 */
	Time t_max = MIN_TIME;
	Voltage u_max = Voltage(std::numeric_limits<Real>::lowest());

	/**
	 * Called whenever the integrator processed a time slice. Checks whether the
	 * current membrane potential is larger than the already seen membrane
	 * potential. If yes, updates the locally stored maximum time and voltage.
	 */
	template <typename State, typename System>
	void record(Time t, const State &s, const System &sys)
	{
		const Voltage u = sys.ode().voltage(s, sys);
		if (u > u_max) {
			u_max = u;
			t_max = t;
		}
	}
};

/**
 * Lowest recursion level of the MultiRecorder class, just does nothing.
 */
template <typename... Recorders>
class MultiRecorder : public NullRecorder {
};

/**
 * Class used to cascade a number of Recorders. Use the make_multi_recorder()
 * method to conveniently construct a Recorder consisting of multiple recorders.
 */
template <typename Recorder, typename... Recorders>
class MultiRecorder<Recorder, Recorders...> : MultiRecorder<Recorders...> {
private:
	Recorder &recorder;

public:
	MultiRecorder(Recorder &recorder, Recorders &... rs)
	    : MultiRecorder<Recorders...>(rs...), recorder(recorder)
	{
	}

	template <typename State, typename System>
	void record(Time t, const State &s, const System &sys)
	{
		recorder.record(t, s, sys);
		MultiRecorder<Recorders...>::record(t, s, sys);
	}
};

/**
 * Helper function which can be used to conveniently create a MultiRecorder
 * instance. Simply pass references to all recorders that should be used to the
 * method and store the result in an auto variable.
 */
template <typename... Recorders>
MultiRecorder<Recorders...> make_multi_recorder(Recorders &... rs)
{
	return MultiRecorder<Recorders...>(rs...);
}
}

#endif /* CINDER_ODE_RECORDER_HPP */

