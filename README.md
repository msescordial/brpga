# route-network-genetic-algo

The Matlab script files above generate the "optimal" _s_ routes of a road network using Genetic Algorithm, where _s_ is predefined no. of routes. To generate candidate routes for initializing the route network, k-Shortest Path and Dijkstra's algorithms are used, in which the source code came from https://www.mathworks.com/matlabcentral/fileexchange/32513-k-shortest-path-yen-s-algorithm and the file is kSP.m.

The main file is mainGAwithTerminals.m.

Two road networks are used:
1. A benchmark network which is the network by Mandl, in which the details are in the following file: network_Mandl.m
2. A real-life road network of a city in India, from the work of Suman and Bolia (2019), in which the details are in the following file: network_india.m

The following user-defined parameters can be changed:
1. _s_ - the number of routes (under the files: network_Mandl.m or network_india.m)
2. _the road network_, either network_mandl() or network_india() (under the file: mainGAwithTerminals.m., line 10)
3. _TerminalNodes_ - the nodes which serve as terminals (under the files: network_Mandl.m or network_india.m)
4. _transfer_time_ - the time required to make one transfer (under the files: network_Mandl.m or network_india.m)
5. _maxiter_ - the maximum iterations (or the number of generations) for genetic algorithm (under the file: mainGAwithTerminals.m., line 6)
6. _population_size_ - the population size for each generation for genetic algorithm (under the file: mainGAwithTerminals.m., line 7)
7. _P_ce_ - the inter-crossover probability (under the file: GeneticAlgo.m, line 133)
8. _P_ca_ - the intra-crossover probability (under the file: GeneticAlgo.m, line 192)
9. _P_m_ - the mutation probability (under the file: GeneticAlgo.m, line 226)

Below are the descriptions for each file:
1. BusRoute.m - finalizes the substring representation of a route
2. GeneticAlgo.m - main code for genetic algorithm
3. GeneticAlgoforSAGA. m - main code for genetic algorithm under simulated annealing 
4. InitialPopulationforGA.m - generates the initial population for genetic algorithm
5. ObjFuncVal.m - code for the objective function value
6. StringtoRouteSet.m - converts the string representation to a set of routes
7. TotalTime.m - outputs the total time matrix of the solution
8. case1feasible.m - converts the infeasible route to a feasible one using Case 1 of the method from the work of (Buba and Lee, 2016)
9. case2feasible.m - converts the infeasible route to a feasible one using Case 2 of the method from the work of (Buba and Lee, 2016)
10. checkAllNodes.m - checks if all nodes are in the solution; outputs 1 if true and 0 if otherwise
11. common_nodes.m - outputs the common nodes of two routes
12. displayBusRoute.m - displays a route
13. displayPopulationSolution.m - displays the string representation of the population of solutions for genetic algorithm
14. displayRouteSet.m - displays a set of routes
15. displaySolution.m - displays the string representation of a solution
16. generateBusRoutes.m - generates the candidate routes
17. generateBusRoutesWithTerminals.m - generates the candidate routes where the terminals are specified
18. generateInitialNetwork.m - generates an initial solution
19. intra_crossover_main.m - main code for intra-crossover
20. intra_crossover_operation.m - subcode for intra-crossover operation
21. kSP.m - generates the shortest paths using k-Shortest Path and Dijkstra's algorithms
22. mainGA.m - main code for genetic algorithm where the terminals are not specified
23. mainGAWithTerminals.m - main code for genetic algorithm where the terminals are specified
24. mainSAGA.m - main code for genetic algorithm and simulated annealing where the terminals are not specified
25. mutation_op.m - main code for mutation operation
26. network_india.m - contains the details of the network of a city in India, from the work of Suman and Bolia (2019)
27. network_mandl.m - contains the details of Mandl's network
28. node_position.m - outputs the position of a node in a route
29. repair_infeasibility.m - main code for repairing an infeasible solution
30. testing.m - file used for testing some codes
31. tijCase1.m - computes the travel time between two nodes in a network where no transfer is needed
32. tijCase2.m - computes the travel time between two nodes in a network where one transfer is needed
33. tijCase3.m - computes the travel time between two nodes in a network where two transfers are needed
