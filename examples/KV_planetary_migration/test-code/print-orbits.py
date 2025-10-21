# Import the rebound module
import rebound
import os
import numpy as np



sa = rebound.SimulationArchive('simulationarchive.bin')
#sa = rebound.SimulationArchive('planet-ic.bin')
outfile = open("orbits.txt","w")
outfile.write("#id, time (years), a, e, i, node, arg peri, mean anomaly (all radians)\n")
for i, sim in enumerate(sa):
    n = sim.N
    com = sim.calculate_com()
    time = sim.t/(2.*np.pi)
    for j in range(n-1):
        p = sim.particles[j+1]
        o = p.calculate_orbit(com)
        outfile.write("%d\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n" % (p.hash.value,time,o.a,o.e,o.inc, o.Omega,o.omega,o.M))
outfile.close()

quit()

sa = rebound.SimulationArchive('high-res-end-archive.bin')
#sa = rebound.SimulationArchive('planet-ic.bin')
outfile = open("high-res-end-orbits.txt","w")
outfile.write("#id, time (years), a, e, i, node, arg peri, mean anomaly (all radians)\n")
for i, sim in enumerate(sa):
    n = sim.N
    com = sim.calculate_com()
    time = sim.t/(2.*np.pi)
    for j in range(n-1):
        p = sim.particles[j+1]
        o = p.calculate_orbit(com)
        outfile.write("%d\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n" % (j+1,time,o.a,o.e,o.inc, o.Omega,o.omega,o.M))
outfile.close()


