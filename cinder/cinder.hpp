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
 * @file cinder.hpp
 *
 * Cinder main header. Includes all models and utility classes. Include
 * individual headers for faster compilation time.
 *
 * @author Andreas Stöckel
 */

#pragma once

#ifndef CINDER_HPP
#define CINDER_HPP

#include <cinder/common/fast_math.hpp>
#include <cinder/common/spikes.hpp>
#include <cinder/common/time.hpp>
#include <cinder/common/types.hpp>
#include <cinder/common/vector.hpp>
#include <cinder/models/current_source.hpp>
#include <cinder/models/neuron.hpp>
#include <cinder/models/synapse.hpp>
#include <cinder/models/neurons/adex.hpp>
#include <cinder/models/neurons/hodgkin_huxley.hpp>
#include <cinder/models/neurons/izhikevich.hpp>
#include <cinder/models/neurons/lif.hpp>
#include <cinder/models/neurons/mat.hpp>
#include <cinder/models/synapses/cond_alpha.hpp>
#include <cinder/models/synapses/cond_exp.hpp>
#include <cinder/models/synapses/cur_alpha.hpp>
#include <cinder/models/synapses/cur_exp.hpp>
#include <cinder/models/synapses/delta.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/integrator/euler.hpp>
#include <cinder/integrator/midpoint.hpp>
#include <cinder/integrator/runge_kutta.hpp>
#include <cinder/ode/controller.hpp>
#include <cinder/ode/ode.hpp>
#include <cinder/ode/recorder.hpp>
#include <cinder/ode/solver.hpp>

#endif /* CINDER_HPP */
