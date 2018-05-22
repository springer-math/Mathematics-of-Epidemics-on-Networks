# Tests for EoN package:  #

## Functional Tests   

Normal functional tests are written using nosetests package. To install: 

	pip install nose 

Nose allow many different way of running tests. For e.g. test names specified may be file or module names, and may optionally indicate the test case to run by separating the module or file name from the test case name with a colon. 

Examples:

	nosetests -v tests\test_sim_sweep_parameters.py
	nosetests -v EoN.tests.test_sim_sweep_parameters
	nosetests -v EoN.tests.test_sim_sweep_parameters:TestSimSweepParameters.test_Gillespie_SIS_type
	
Example test module execution with verbose logs: 

	C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\tests>nosetests -v test_sim_sweep_parameters.py
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_Gillespie_SIR_sweep_gamma ... ok
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_Gillespie_SIR_sweep_tau ... ok
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_Gillespie_SIS_sweep_gamma ... ERROR
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_Gillespie_SIS_sweep_tau ... ERROR
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_Gillespie_SIS_type ... ok
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_basic_discrete_SIR_sweep_p ... ok
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_basic_discrete_SIS_sweep_p ... ok
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_fast_SIR_sweep_gamma ... ok
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_fast_SIR_sweep_tau ... ok
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_fast_SIS_sweep_gamma ... ok
	EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_fast_SIS_sweep_tau ... ok
	
	======================================================================
	ERROR: EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_Gillespie_SIS_sweep_gamma
	----------------------------------------------------------------------
	Traceback (most recent call last):
	  File "c:\users\tting\appdata\local\programs\python\python36\lib\site-packages\nose\case.py", line 198, in runTest
	    self.test(*self.arg)
	  File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\tests\test_sim_sweep_parameters.py", line 41, in test_Gillespie_SIS_sweep_gamma
	    [EoN.Gillespie_SIS(G, 1.0, i, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]
	  File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\tests\test_sim_sweep_parameters.py", line 41, in <listcomp>
	    [EoN.Gillespie_SIS(G, 1.0, i, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]
	  File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\simulation.py", line 3039, in Gillespie_SIS
	    node_history = _transform_to_node_history_(infection_times, recovery_times, tmin, SIR = False)
	  File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\simulation.py", line 136, in _transform_to_node_history_
	    Rtimes = recovery_times[node]
	KeyError: (24, 24)
	
	======================================================================
	ERROR: EoN.tests.test_sim_sweep_parameters.TestSimSweepParameters.test_Gillespie_SIS_sweep_tau
	----------------------------------------------------------------------
	Traceback (most recent call last):
	  File "c:\users\tting\appdata\local\programs\python\python36\lib\site-packages\nose\case.py", line 198, in runTest
	    self.test(*self.arg)
	  File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\tests\test_sim_sweep_parameters.py", line 53, in test_Gillespie_SIS_sweep_tau
	    [EoN.Gillespie_SIS(G, i, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]
	  File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\tests\test_sim_sweep_parameters.py", line 53, in <listcomp>
	    [EoN.Gillespie_SIS(G, i, 1.0, initial_infecteds=initial_infections, return_full_data=True, tmax=10) for i in np.arange(0.0, 1.0, 0.1)]
	  File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\simulation.py", line 3039, in Gillespie_SIS
	    node_history = _transform_to_node_history_(infection_times, recovery_times, tmin, SIR = False)
	  File "C:\GitRepos\tinghf_Math_Epidemics_Networks\EoN\simulation.py", line 136, in _transform_to_node_history_
	    Rtimes = recovery_times[node]
	KeyError: (24, 24)
	
	----------------------------------------------------------------------
	Ran 11 tests in 2.435s
	
	FAILED (errors=2)

For more details on usage, please see [http://nose.readthedocs.io/en/latest/usage.html](http://nose.readthedocs.io/en/latest/usage.html)
 

## Performance Tests  

Perf tests are created using the Python perf modules. To install: 

	pip install perf

It features: 

- Simple API to run reliable benchmarks
- Automatically calibrate a benchmark for a time budget.
- Spawn multiple worker processes.
- Compute the mean and standard deviation.
- Detect if a benchmark result seems unstable.
- JSON format to store benchmark results.
- Support multiple units: seconds, bytes and integer. 

To run a benchmark use the perf timeit command, with result written into perf_run.json, as follows: 

	python tests\perf_run.py -o perf_run.json

An example to write your own benchmark script bench.py, is as follows: 

	#!/usr/bin/env python3
	import perf
	
	runner = perf.Runner()
	runner.timeit(name="sort a sorted list",
	              stmt="sorted(s, key=f)",
	              setup="f = lambda x: x; s = list(range(1000))")

Then you could take perf stats on this script simply by: 

	python bench.py -o bench.json
	.....................
	sort a sorted list: Mean +- std dev: 93.4 us +- 7.4 us

The benchmark file generated for each execution include the environment context(so you could do a apple-2-apple comparison of benchmarks later on). It is at the tail end of the json file generated and can also be displayed by: 

	python -m perf metadata bench.json
	Metadata:
	- cpu_count: 4
	- hostname: NA00327L
	- loops: 2048
	- name: sort a sorted list
	- perf_version: 1.5.1
	- platform: Windows-10-10.0.14393-SP0
	- python_executable: C:\Users\tting\AppData\Local\Programs\Python\Python36\python.exe
	- python_implementation: cpython
	- python_version: 3.6.4 (64-bit) revision d48eceb
	- timeit_setup: 'f = lambda x: x; s = list(range(1000))'
	- timeit_stmt: 'sorted(s, key=f)'
	- timer: QueryPerformanceCounter(), resolution: 395 ns
	- unit: second

To analyze benchmark results use the perf stats command:

	python -m perf stats bench.json
	Total duration: 16.6 sec
	Start date: 2018-05-18 16:23:58
	End date: 2018-05-18 16:24:18
	Raw value minimum: 178 ms
	Raw value maximum: 234 ms
	
	Number of calibration run: 1
	Number of run with values: 20
	Total number of run: 21
	
	Number of warmup per run: 1
	Number of value per run: 3
	Loop iterations per value: 2048
	Total number of values: 60
	
	Minimum:         87.0 us
	Median +- MAD:   90.8 us +- 1.9 us
	Mean +- std dev: 93.4 us +- 7.4 us
	Maximum:         114 us
	
	  0th percentile: 87.0 us (-7% of the mean) -- minimum
	  5th percentile: 87.7 us (-6% of the mean)
	 25th percentile: 89.2 us (-4% of the mean) -- Q1
	 50th percentile: 90.8 us (-3% of the mean) -- median
	 75th percentile: 93.4 us (-0% of the mean) -- Q3
	 95th percentile: 113 us (+22% of the mean)
	100th percentile: 114 us (+23% of the mean) -- maximum
	
	Number of outlier (out of 83.0 us..99.6 us): 7

Or you could render an histogram in text mode as follows: 

	python -m perf hist bench.json
	86.5 us:  3 ################
	87.6 us:  5 ##########################
	88.7 us: 15 ###############################################################################
	89.8 us:  7 #####################################
	90.9 us: 11 ##########################################################
	92.0 us:  3 ################
	93.1 us:  2 ###########
	94.2 us:  2 ###########
	95.3 us:  1 #####
	96.4 us:  1 #####
	97.5 us:  2 ###########
	98.5 us:  1 #####
	99.6 us:  1 #####
	 101 us:  0 |
	 102 us:  0 |
	 103 us:  0 |
	 104 us:  0 |
	 105 us:  0 |
	 106 us:  0 |
	 107 us:  0 |
	 108 us:  0 |
	 109 us:  0 |
	 111 us:  0 |
	 112 us:  1 #####
	 113 us:  2 ###########
	 114 us:  3 ################

Use compare_to to compare between benchmarks to see if there's no significant difference:
 
	python -m perf compare_to bench.json bench2.json
	Benchmark hidden because not significant (1): sort a sorted list

or not 

	python -m perf compare_to bench.json bench2.json
	+-----------+---------+------------------------------+
	| Benchmark | bench   | bench2                       |
	+===========+=========+==============================+
	| timeit    | 4.70 us | 4.22 us: 1.11x faster (-10%) |
	+-----------+---------+------------------------------+


For more details on usage, please see [https://github.com/vstinner/perf](https://github.com/vstinner/perf)
 
