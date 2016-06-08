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

#include <cmath>
#include <fstream>
#include <functional>

#include "gtest/gtest.h"

#include <cinder/models/neurons/mat.hpp>
#include <cinder/models/current_source.hpp>
#include <cinder/integrator/dormand_prince.hpp>
#include <cinder/ode/solver.hpp>
#include <cinder/ode/recorder.hpp>

namespace cinder {
TEST(mat, constant_current)
{
	// Expected spike times as calculated with the MAT demo program with
	// increased time step (0.01ms) and fixed differential equation (R/tau)
	// instead of R. See:
	// http://www.ton.scphys.kyoto-u.ac.jp/~shino/toolbox/kspred/pred.html
	const std::vector<Time> expected_spikes = {
	    2.379998_ms,   8.530100_ms,   16.610285_ms,  25.680492_ms,
	    35.290085_ms,  45.328400_ms,  55.786644_ms,  66.665833_ms,
	    77.988251_ms,  89.730759_ms,  101.913361_ms, 114.506050_ms,
	    127.518829_ms, 140.921829_ms, 154.694260_ms, 168.826492_ms,
	    183.288544_ms, 198.060425_ms, 213.112152_ms, 228.413742_ms,
	    243.945206_ms, 259.672150_ms, 275.597687_ms, 291.673370_ms,
	    307.879181_ms, 324.205109_ms, 340.621124_ms, 357.127228_ms,
	    373.693390_ms, 390.329620_ms, 407.005890_ms, 423.732208_ms,
	    440.488556_ms, 457.274933_ms, 474.091339_ms, 490.927765_ms,
	    507.774200_ms, 524.640625_ms, 541.517090_ms, 558.403564_ms,
	    575.300049_ms, 592.206543_ms, 609.113037_ms, 626.029541_ms,
	    642.946045_ms, 659.862549_ms, 676.789062_ms, 693.715576_ms,
	    710.642090_ms, 727.568604_ms, 744.495117_ms, 761.431641_ms,
	    778.368164_ms, 795.304688_ms, 812.241211_ms, 829.177734_ms,
	    846.114258_ms, 863.050781_ms, 879.987305_ms, 896.923828_ms,
	    913.860352_ms, 930.796875_ms, 947.733398_ms, 964.669922_ms,
	    981.606445_ms, 998.542969_ms};

	DormandPrinceIntegrator integrator;
	NullRecorder recorder;
	NullController controller;
	std::vector<Time> spikes;

	// Create a linear integrate and fire neuron with a current based synapse
	auto current_source = ConstantCurrentSource(5_nA);
	auto neuron = make_neuron<MAT2>(current_source,
	                                [&spikes](Time t) { spikes.push_back(t); });

	// Solve the equation
	make_solver(neuron, integrator, recorder, controller).solve(1_s);

	// Make sure the two results are within the range of 0.3ms
	ASSERT_EQ(expected_spikes.size(), spikes.size());
	for (size_t i = 0; i < spikes.size(); i++) {
		EXPECT_GT(0.3_ms, std::abs(expected_spikes[i] - spikes[i]));
	}
}
}
