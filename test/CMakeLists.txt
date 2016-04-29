#  Cinder -- C++ Single Spiking Neuron Simulator
#  Copyright (C) 2015, 2016  Andreas Stöckel, Christoph Jenzen
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# CMake version and project metadata
cmake_minimum_required(VERSION 3.2)
project(cinder)

include_directories(
	PRIVATE
		${GTEST_INCLUDE_DIRS}
)

add_executable(test_cinder
	common/test_vector
	models/test_current_source
	models/neurons/test_adex
	models/neurons/test_hodgkin_huxley
	models/neurons/test_izhikevich
	models/neurons/test_lif
	models/synapses/test_cur_alpha
	models/synapses/test_cur_exp
	models/synapses/test_delta
	ode/test_solver
	ode/test_controller
)
target_link_libraries(test_cinder
	cinder
	-pthread
	${GTEST_LIBRARIES}
)
add_test(test_cinder test_cinder)

