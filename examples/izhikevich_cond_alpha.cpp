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

