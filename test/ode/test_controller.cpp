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

#include "gtest/gtest.h"

#include <cinder/ode/controller.hpp>

namespace cinder {

namespace {
template<typename...T>
static void test_multi_controller(ControllerResult expected, const T&... ts) {
	auto mc = make_multi_controller(ts...);
	EXPECT_EQ(expected, mc.control(0_s, 0, 0));
}
}

TEST(controller, multi_controller)
{
	ConstantController<ControllerResult::CONTINUE> cc;
	ConstantController<ControllerResult::MAY_CONTINUE> cmc;
	ConstantController<ControllerResult::ABORT> ca;

	test_multi_controller(ControllerResult::CONTINUE, cc);
	test_multi_controller(ControllerResult::MAY_CONTINUE, cmc);
	test_multi_controller(ControllerResult::ABORT, ca);

	test_multi_controller(ControllerResult::CONTINUE, cc, cc);
	test_multi_controller(ControllerResult::CONTINUE, cc, cmc);
	test_multi_controller(ControllerResult::ABORT, cc, ca);

	test_multi_controller(ControllerResult::CONTINUE, cmc, cc);
	test_multi_controller(ControllerResult::MAY_CONTINUE, cmc, cmc);
	test_multi_controller(ControllerResult::ABORT, cmc, ca);

	test_multi_controller(ControllerResult::ABORT, ca, cc);
	test_multi_controller(ControllerResult::ABORT, ca, cmc);
	test_multi_controller(ControllerResult::ABORT, ca, ca);

	test_multi_controller(ControllerResult::CONTINUE, cc, cc, cc);
	test_multi_controller(ControllerResult::CONTINUE, cc, cc, cmc);
	test_multi_controller(ControllerResult::ABORT, cc, cc, ca);
	test_multi_controller(ControllerResult::CONTINUE, cc, cmc, cc);
	test_multi_controller(ControllerResult::CONTINUE, cc, cmc, cmc);
	test_multi_controller(ControllerResult::ABORT, cc, cmc, ca);
	test_multi_controller(ControllerResult::ABORT, cc, ca, cc);
	test_multi_controller(ControllerResult::ABORT, cc, ca, cmc);
	test_multi_controller(ControllerResult::ABORT, cc, ca, ca);

	test_multi_controller(ControllerResult::CONTINUE, cmc, cc, cc);
	test_multi_controller(ControllerResult::CONTINUE, cmc, cc, cmc);
	test_multi_controller(ControllerResult::ABORT, cmc, cc, ca);
	test_multi_controller(ControllerResult::CONTINUE, cmc, cmc, cc);
	test_multi_controller(ControllerResult::MAY_CONTINUE, cmc, cmc, cmc);
	test_multi_controller(ControllerResult::ABORT, cmc, cmc, ca);
	test_multi_controller(ControllerResult::ABORT, cmc, ca, cc);
	test_multi_controller(ControllerResult::ABORT, cmc, ca, cmc);
	test_multi_controller(ControllerResult::ABORT, cmc, ca, ca);

	test_multi_controller(ControllerResult::ABORT, ca, cc, cc);
	test_multi_controller(ControllerResult::ABORT, ca, cc, cmc);
	test_multi_controller(ControllerResult::ABORT, ca, cc, ca);
	test_multi_controller(ControllerResult::ABORT, ca, cmc, cc);
	test_multi_controller(ControllerResult::ABORT, ca, cmc, cmc);
	test_multi_controller(ControllerResult::ABORT, ca, cmc, ca);
	test_multi_controller(ControllerResult::ABORT, ca, ca, cc);
	test_multi_controller(ControllerResult::ABORT, ca, ca, cmc);
	test_multi_controller(ControllerResult::ABORT, ca, ca, ca);
}
}

