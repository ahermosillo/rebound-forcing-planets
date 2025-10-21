explanation for input files

units for input files: au, yr, deg
units in code: au, 2pi*yr, rad

====== integration-parameters.txt =====
integration timestep, total integration time, output timestep, trial sim (0 or 1), number of planets

trial sim 0 = test particles will be added (see tp-parameters.txt)
trial sim 1 = no test particles added, just testing planet behavior

in the code, the time is converted to yr 2 pi units

====== planet-ic.txt =======
mass, radius, semimajor axis, eccentricity, inclination, argument of pericenter, longitude of ascending node, true anomaly 

repeats again for each planet.

====== planet-migration-parameters.txt =========
tau_a, tau_e, tau_i, tau_w, tau_node, final_a, final_e, final_i, final_w, final_node

repeat again for each planet. 

final_a - initial_a = magnitude of change for desired evolution

====== tp-parameters.txt ==========
number of particles, min semimajor axis, max semimajor axis, delta e, delta i, randomized seed

delta_e and delta_i are assuming a raleigh distribution. 
