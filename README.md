# multiplexing-noise
Open source code for paper "The synchronizing role of multiplexing noise in networks of Kuramoto phase oscillators"

This code is usable for simulating a two layer network of Kuramoto phase oscillators, both with and without noise. The interlayer coupling is multiplex, meaning each oscillator in one layer is connected to a single oscillator in the other layer. The coupling within the layer is non-local, meaning that each oscillator is connected to the R nearest oscillators on each side.  

The simulations are performed with C++. The integrator in the main cpp files require two inputs, namely the stepper and the system. The stepper files are the integration methods (for example Runge-Kutta 4th order) and the system files are the differential equations themselves. The observer files are used to observe the output of the integrator, which can be output to a data file. The noise files are for generating colored noise. The 'det' files are for the deterministic system and use the standard Runge-Kutta method. The 'stoch' files are for simulating the system with white noise and use the stochastic Runge-Kutta method. Fractional Brownian motion is integrated using the Euler-Maruyama method. The tplot files are for outputting the oscillator phases and phase frequencies. The state files output the averaged order parameters after a transient time has elapsed. Plotting is performed with python. 

If the code is helpful for research then please consider citing the paper "The synchronizing role of multiplexing noise in networks of Kuramoto phase oscillators".  
