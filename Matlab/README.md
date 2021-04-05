This orbit propagator began development in AA 279A and was additionally supplemented.

Matlab version: 2017b.

Top-level file descriptions:
- run_single_sim
  - Allows user to run a simulation then plots + displays sim results.
- Propagator.slx
  - The guts of the simulation.
  - IMPORTANT NOTE: Many calculations assume e = 0 = constant. Probably should be updated.
Subdirectory descriptions:
- functions_and_helper_scripts
  - Contains matlab functions and scripts that are called by the propagator and top-level scripts.
  - Some internal functions are stored within Propagator.slx