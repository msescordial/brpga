clc
clear

tic

% Inputs for GA
max_no_of_generations = 20;         % maximum no. of generations
population_size = 100;              % must be divisible by 4
P_ce = 0.5;             % inter-crossover probability
P_ca = 0.5;             % intra-crossover probability
P_m = 0.05;             % mutation probability

% Constraints
min_route_length = 3;

% Choose Network
[DistanceMatrix,TimeMatrix,TravelDemandMatrix,TerminalNodes,k_kSP,s,transfer_time]=network_mandl();
n = size(DistanceMatrix,1);         % n = no. of nodes

% Genetic Algorithm    
[S] = EliteGeneticAlgorithm(max_no_of_generations, population_size, P_ce, P_ca, P_m, min_route_length, ...
            DistanceMatrix, TimeMatrix, TravelDemandMatrix, TerminalNodes, k_kSP, s, transfer_time, n);
[Sr] = stringToRoutes(S,s,n); 

% Optimal Route Set and its Objective Function Value
fprintf('\n\n Optimal Route Set S*: \n '); 
displayRouteSet(S,s,n);
for a=1:s
    fprintf(' Route %d:', a); 
    br = functionRoute(Sr{a,1});
    displayRoute(br);
end
SolutionTimeMatrix = getRouteSetTimeMatrix(Sr,s,TimeMatrix, transfer_time);
%fprintf("\n Route Set Time Matrix \n"); disp(SolutionTimeMatrix);
E = getObjectiveFunctionValue(Sr,TravelDemandMatrix,DistanceMatrix,SolutionTimeMatrix,n);   
fprintf('\n Objective Function Value E*: %f \n\n', E);

toc