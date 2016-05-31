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

#include <cinder/common/spikes.hpp>

namespace cinder {
TEST(spikes, poisson_spike_train)
{
	for (float lambda = 0.0; lambda < 20.0; lambda += 0.1) {
		const int seed = 173812;
		const std::vector<Time> spike_train =
			poisson_spike_train(0_s, 2000_s, lambda, seed);
		EXPECT_NEAR(lambda, float(spike_train.size()) / 2000.0f, 0.1);
	}
}
}

