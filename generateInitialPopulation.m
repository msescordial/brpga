function init_pop_matrix = generateInitialPopulation(population_size, S0r, s, n, BusRouteID, TotalNoOfRoutes, ...
        DistanceMatrix, TimeMatrix, TravelDemandMatrix, waiting_time, transfer_time)

    % Initial Population for GA

    init_pop_matrix = cell(population_size,3);
    
    % Creating the Initial Population 
    % 1st col: no., 2nd col: route nos., 3rd col: Obj Func Value, 4th col: Fitness Value
    % First Row (S0)
    init_pop_matrix{1,1} = 1;
    init_pop_matrix{1,2} = S0r;             %init_pop_matrix{1,2} = transpose(S0);
    S0_SolutionTimeMatrix = TotalTime(S0r,s,TimeMatrix, waiting_time, transfer_time);
    ofv0 = ObjFuncVal(S0r,TravelDemandMatrix,TimeMatrix,S0_SolutionTimeMatrix,n);
    init_pop_matrix{1,3} = ofv0;   init_pop_matrix{1,4} = ofv0;   

    % Succeeding Rows
    for g = 1: population_size-1
        init_pop_matrix{g+1,1} = g+1;
        [S1,S1r] = generateInitialNetwork(DistanceMatrix, BusRouteID, TotalNoOfRoutes,s);
        init_pop_matrix{g+1,2} = S1r;
        S1_SolutionTimeMatrix = TotalTime(S1r,s,TimeMatrix, waiting_time, transfer_time);
        ofv1 = ObjFuncVal(S1r,TravelDemandMatrix,TimeMatrix,S1_SolutionTimeMatrix,n);
        init_pop_matrix{g+1,3} = ofv1; init_pop_matrix{g+1,4} = ofv1;  
    end
    % disp("Initial Population Matrix"); disp(init_pop_matrix);

end
