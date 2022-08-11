# Heat_Microstructure
This repo simulate the heat transfer and microstructure evolution during laser powder bed fusion, which is an additive manufacturing process.
This repo includes 3 independent sub-repos: CA_Heat, compare_analytical_FVM, heat_adaptive_moltenTrack.

## CA_Heat
This is a c-based, process-oriendted, serial program that simulates the heat transfer and microstructure in LPBF. This program is used for illustrating the basic physical mechanism used for simulation, such as heat transfer simulation using simplest analytical method, along with the dendrite growth simultion using cellular automata algorithm. This program did not apply any complex and adaptive numerical methods. Due to limitation of of serial program, this program can only deal with the simulation with very large time step and mesh size, and very short total time.

## compare_analytical_FVM
This is a framework used for benchmarking performance and comparing consistency (simulation results vs experiment results) for different numerical methods. 

This cpp-based program also allow user to define or create complex scanning strategy of laser scanning pattern, such as spiral pattern, zig-zag pattern, etc. User can also choose to overwrite the laser position and laser speed function to customized their own scanning pattern with any path shape and velocities. 

In the tempClass.h, implemented mainstream numerical methods includes analytical method, semi-analytical method, finite volume method with isolated boundary condition, finite volume method with radiative and convective boundary condition, semi-analytical method with mirroring boundary condition. 

This program output a serices of .dat file. Each .dat file is the temperature field of 3d field at a specific time. This program also have a post-process program written in python to extract, plot, and compare the temperature of several points predicted using different numerical methods.

## heat_adaptive_molten_pool
This is the final framework that use semi-analytical method with mirroring boundary condition for heat transfer method, along with cellular automata for microstructure prediction method. 

This program implemented a molten-pool-track algorithm that dynamically track the molten pool points during the previous time step, and quickly check and update the temperature within the motlen pool at current time step in a "envoling" way. This algorithm allow us to track the temperature from the laser spot center (which is highest temp) downwards to any customized temerature contours. In this way, the heat transfer in additive manufacturing, which usually has a relatively high spatial locality, could be accelerated without significantly losing accuracy.



