# rebound-forcing-planets
This repo contains modified REBOUND code which allows the user to evolve particles' orbital elements given arbitrary functions. 

The original REBOUND code can be found [here](https://github.com/hannorein/rebound) and documentation can be found [here](https://rebound.readthedocs.io/en/latest/). 

My modifications can be found in srcDev and only for the mercurius and whfast integrators. The modifications were made for REBOUND version 3.17.2 and then updated to 3.19.2 to fix a bug in the mercurius integrator. 

Full details can be found in ApJ paper by Hermosillo Ruiz et. al (2025) titled “Forcing Planets to Evolve: Interactions Between Uranus and Neptune at Late Stages of Dynamical Evolution.” [see paper](https://ui.adsabs.harvard.edu/abs/2025ApJ...987..125H/abstract)

The simulations, simulation output, and plotting scripts used to create each figure in the manuscript can be found in HermosilloRuiz2025_examples.
