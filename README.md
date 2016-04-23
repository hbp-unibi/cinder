# Cinder – A hyper-fast single spiking neuron simulator

[![Build Status](https://travis-ci.org/hbp-sanncs/cinder.svg?branch=master)](https://travis-ci.org/hbp-sanncs/cinder)

## About

Spiking neuron models are generally given as a set of differential equations
which usually cannot be solved in closed form. Unfortunately, this implies that
analysing the behaviour of these models is – in most cases – restricted to
numerical simulation of the neuron dynamics.

_Cinder_ is a modular library which implements numerical simulation of single
neurons, allowing to analyse the properties of the neuron models. In contrast
to network simulators _Cinder_ gives full control over the simulation process.
It allows to record any neuron state variable, to use any number of current
sources (such as synapses) and to use multiple different integrators.

## Installing and Requirements

_Cinder_ requires a C++14 compatible compiler such as `GCC 5` or `clang 3.6`,
as well as `cmake 3.2`. You can build and install _Cinder_ using the following
commands:

```bash
git clone --depth 1 https://github.com/hbp-sanncs/cinder && cd cinder
mkdir build && cd build
cmake .. && make && make test
sudo make install
```

## Example Code

The following code simulates an Izhikevich neuron with two conductance based
synapses – one excitatory synapse and one inhibitory synapse. The excitatory
synapse receives spikes at about 100 and 500 ms, the second synapse at about
400 ms into the simulation.


```c++
#include <iostream>

#include <cinder/cinder.hpp>

using namespace cinder;

int main()
{
	// Use the DormandPrince integrator with default target error
	DormandPrinceIntegrator integrator;

	// Record the neuron state as CSV to cout
	CSVRecorder recorder(std::cout);

	// Use the NeuronController class to automatically abort the simulation
	// once the neuron has settled to its resting state
	NeuronController controller;

	// Assemble current source (here: two synapses)
	auto current_source = make_current_source(
	    CondAlpha(20_nS, 10_ms, 30_mV, {95_ms, 100_ms, 102_ms, 110_ms, 499_ms,
	                                     505_ms, 510_ms, 520_ms, 530_ms}),
	    CondAlpha(50_nS, 30_ms, -80_mV, {390_ms, 400_ms, 420_ms, 450_ms}));

	// Assemble the neuron (an Izhikevich neuron with default parameters)
	auto neuron = make_neuron<Izhikevich>(current_source, [] (Time) {
		// Lambda called whenever the neuron spikes
	});

	// Simulate the resulting ODE -- the recorder will record the results
	make_solver(neuron, integrator, recorder, controller).solve();
}
```

Plotting the resulting membrane potential and conductivity traces results in
the following image:

![Izhikevich Cond Alpha Results](https://raw.githubusercontent.com/hbp-sanncs/cinder/master/docs/izhikevich_cond_alpha.png "Result of the Izhikevich-Cond-Alpha example program")


## Features

_Cinder_ is implemented as a C++ almost-header-only template library. This
allows the compiler to construct a customised code-path for exactly the
simulation that is being performed, while keeping the actual model
implementations extremely concise. With enabled optimisation, a single neuron
simulation is converted into an inline blob of vectorised assembly without any
external function calls. This way, _Cinder_ is magnitudes faster than full
network simulators simulating a single neuron.

Conceptually, each single neuron simulation is represented by an autonomous
ordinary differential equation (ODE) which is automatically assembled at compile
time using variadic templates. This single ODE representing the entire neuron,
including its membrane and current sources, is then solved using a numerical
differential equation integrator.

### Implemented neuron models

* _LIF_: Linear Integrate and Fire
* _AdEx_: Exponential Integrate and Fire Neuron with adaptive threshold
* _Izhikevich_: Quadratic Integrate and Fire Neuron

New neuron models can be easily implemented by simply defining the parameter and
state vector and the corresponding differential equation. More neuron models
will be added in the future.

### Implemented synapse models

* _CondAlpha_: Conductance-based synapse with alpha-function shaped decay
* _CondExp_: Conductance-based synapse with exponential-function shaped decay
* _CurAlpha_: Current-based synapse with alpha-function shaped decay
* _CurExp_: Current-based synapse with exponential-function shaped decay
* _Delta_: Dirac-Delta current pulse causing a defined voltage step

Furthermore, arbitrary current sources can be used to inject currents into the
neuron. Pre-defined current sources include:

* _ConstantCurrentSource_: Current source which injects a constant current into
the neuron
* _StepCurrentSource_: Current source which injects a constant current during a
pre-defined time span into the neuron

### Implemented integrators

* _EulerIntegrator_: The most simple differential equation integrator
* _MidpointIntegrator_: Second-order Runge-Kutta integrator
* _RungeKuttaIntegrator_: Fourth-order Runge-Kutta integrator
* _DormandPrinceIntegrator_: Fourth-order Runge-Kutta-based adaptive stepsize integrator

## Authors

This project has been initiated by Andreas Stöckel and Christoph Jenzen in 2016 while working
at Bielefeld University in the [Cognitronics and Sensor Systems Group](http://www.ks.cit-ec.uni-bielefeld.de/) which is
part of the [Human Brain Project, SP 9](https://www.humanbrainproject.eu/neuromorphic-computing-platform).

## License

This project and all its files are licensed under the
[GPL version 3](http://www.gnu.org/licenses/gpl.txt) unless explicitly stated
differently.
