function [S1] = EliteGeneticAlgorithm(max_no_of_generations, population_size, P_ce, P_ca, P_m, min_route_length, ...
            DistanceMatrix, TimeMatrix, TravelDemandMatrix, TerminalNodes, k_ksP, s, transfer_time, n)
    
% Binary (1 if to be performed, or 0 if not)
intercrossover = 1;
intracrossover = 1;

elite_ratio = 0.10;      % Elite Population Proportion

% ----- INITIALIZATION -----
% GENERATE CANDIDATE ROUTES
[BusRouteID, AllPaths, AllCosts, TotalNoOfRoutes] = generateRoutes(DistanceMatrix,k_ksP,TerminalNodes);
fprintf('No. of candidate routes generated is %d\n\n', TotalNoOfRoutes);

% GENERATE INITIAL POPULATION (ROUTE IDs)
initial_pop_matrix = cell(population_size,4);
% 1st col: no., 2nd col: route IDs, 3rd col: Obj Func Value, 4th col: Fitness Value
% Rows
for g = 1: population_size
    initial_pop_matrix{g,1} = g;
    [S1,S1r] = generateInitialRouteSet(DistanceMatrix, BusRouteID, TotalNoOfRoutes, s, min_route_length);
    initial_pop_matrix{g,2} = transpose(S1);
    S1_SolutionTimeMatrix = getRouteSetTimeMatrix(S1r,s,TimeMatrix, transfer_time);
    ofv1 = getObjectiveFunctionValue(S1r,TravelDemandMatrix,TimeMatrix,S1_SolutionTimeMatrix,n);
    initial_pop_matrix{g,3} = ofv1; initial_pop_matrix{g,4} = ofv1;  
end
%disp("Initial Population (Route IDs) Matrix"); disp(initial_pop_matrix);

% INITIAL POPULATION (STRINGS)
% Note: Demarcation Line is at the end of every nth node
string_length = s*n;
% Initial Population: Convert from Route ID to String
initial_pop = zeros(population_size, s);
for h=1:population_size
    initial_pop(h,:) = initial_pop_matrix{h,2};
end
%disp("Initial Population (Route IDs):"); disp(initial_pop);
initial_pop_string = zeros(population_size,string_length);
for h=1:population_size
    y1=1; y2=n;
    for k_ksP=1:s
        route_string = BusRouteID{initial_pop(h,k_ksP),2};
        initial_pop_string(h,y1:y2)=route_string;
        y1=y1+n; y2=y2+n;
    end
end
%disp("Initial Population String:"); displayPopulation(initial_pop_string,n);

k = 0;

% --- (GENETIC ALGORITHM) GA LOOP ---
while (k <= max_no_of_generations)

    if (k == 0)
        current_pop_string = initial_pop_string;
    else
        current_pop_string = next_pop_string;
    end

    % EVALUATE ROUTE SETS
    string_objval = cell(population_size,2);       % 1st col: string, 2nd col: objval
    for w=1:population_size
        string_objval{w,1} = current_pop_string(w,:);
    end
    for w=1:population_size
        set_of_routes = cell(s,1);
        p = 1;
        for z=1:s
            set_of_routes{z,1} = current_pop_string(w,p:p+n-1);
            p = p+n;
        end
        %disp("set of routes"); disp(set_of_routes);
        sor_SolutionTimeMatrix = getRouteSetTimeMatrix(set_of_routes,s,TimeMatrix, transfer_time);
        %disp("Solution Time Matrix"); disp(sor_SolutionTimeMatrix);
        sor_objval = getObjectiveFunctionValue(set_of_routes,TravelDemandMatrix,DistanceMatrix,sor_SolutionTimeMatrix,n);
        string_objval{w,2} = sor_objval;
    end
    % Sort; Rank1 = Highest Fitness Value (Minimum Objective Function Value)
    rx = zeros(1,population_size);
    ry = cell(1,population_size);
    for i=1:population_size
        rx(1,i) = string_objval{i,2};
        ry{1,i} = string_objval{i,1};
    end
    [rxsorted,I] = sort(rx);
    rysorted = ry(I);
    %disp("rxsorted"); disp(rxsorted); disp("rysorted"); disp(rysorted); 
    %best_objfuncval = rxsorted(1,1);        
    sorted_pop = zeros(population_size, s*n);
    for j=1:population_size
        sorted_pop(j,:)=rysorted{1,j};
    end
    % disp("Sorted Population:"); displayPopulation(sorted_pop,n);

    %current_pop_string = sorted_pop;

    % Elite Population
    elite_pop = zeros(population_size*elite_ratio, string_length);
    for j=1:population_size*elite_ratio
        elite_pop(j,:)=sorted_pop(j,:);
    end

    next_pop_string = zeros(population_size,string_length);

    % Pass Elite Population to Next Generation
    next_pop_string(1:population_size*elite_ratio,:) = elite_pop;
    
    % display the best solution for each iteration/generation
    S_k = elite_pop(1,:);
    S_k_routes = stringToRoutes(S_k,s,n);
    SolutionTimeMatrix = getRouteSetTimeMatrix(S_k_routes,s,TimeMatrix, transfer_time);
    E_k = getObjectiveFunctionValue(S_k_routes,TravelDemandMatrix,DistanceMatrix,SolutionTimeMatrix,n);   
    fprintf("\n Objective Function Value of Best Route Set of Generation %d: %f",k,E_k); 
    
    % Other Possible Termination Criteria
    %diff = var(rxsorted); disp(diff);  
    %if (diff < 10^(-10))
    %    break;
    %end
    if (k == max_no_of_generations)
        S1 = S_k;
        %disp("S1"); displayRouteSet(S1,s,n);
        break;
    end

    % CROSSOVER
    % Inter-string
    % 1. Rearrangement of Parent population randomly
    % 2. Choosing of demarcation site randomly
    % 3. Crossover Probability P_ce 
    
    if (intercrossover == 1)
    % 1. crossover_matrix has size pop_size*(1-elite_ratio) but still includes elite_pop for crossover, 
    % note: next_pop_string = elite_pop + crossover_matrix
    crossover_matrix1 = zeros(population_size*(1-elite_ratio),string_length);

    % Random pairing of parent solutions for crossover
    r = randperm(population_size);      
    for m=1:population_size*(1-elite_ratio)
        crossover_matrix1(m,:) = sorted_pop(r(m),:);
    end
    %disp("Crossover Rearrangement"); %disp(crossover_matrix1);
    %displayPopulation(crossover_matrix1,n);
    
    % 2. & 3.
    % Demarcation site is at the end of nth element, 2n, 3n, m*(s-1)
    % Only one demarcation site
    
    % Crossover for each pair
    after_intercross = zeros(population_size*(1-elite_ratio),string_length);
    for p1=1:population_size*(1-elite_ratio)*0.5
        strings_for_crossover = zeros(2,s*n);
        strings_for_crossover(1,:)=crossover_matrix1(2*p1-1,:);
        strings_for_crossover(2,:)=crossover_matrix1(2*p1,:);
        %disp("Crossover Strings"); displayRouteSet(crossover_strings(1,:),s,n); displayRouteSet(crossover_strings(2,:),s,n);
        r1 = rand;       %disp("r"); disp(r);
        if ( r1 > P_ce )         % then crossover 
            new_string = zeros(2,s*n);
            demarc_line = n*randi([1,s-1],1) + 1;
            %disp("Demarcation Line"); disp(demarc_line);
            new_string(1,1:demarc_line-1) = strings_for_crossover(1,1:demarc_line-1);
            new_string(1,demarc_line:s*n) = strings_for_crossover(2,demarc_line: s*n);
            new_string(2,1:demarc_line-1) = strings_for_crossover(2,1:demarc_line-1);
            new_string(2,demarc_line:s*n) = strings_for_crossover(1,demarc_line: s*n);
                    
            %disp("Crossover String"); disp(crossover_strings);
            %disp("New String"); disp(new_string);
            
            after_intercross(2*p1-1,:) = new_string(1,:);
            after_intercross(2*p1,:) = new_string(2,:);            
        else
            after_intercross(2*p1-1,:) = strings_for_crossover(1,:);
            after_intercross(2*p1,:) = strings_for_crossover(2,:);
        end
        %disp("After Inter-Crossover");
        %displayRouteSet(after_intercross(2*s1-1,:),s,n);
        %displayRouteSet(after_intercross(2*s1,:),s,n);
    end

    % checking if all nodes are in the solution set
    for z=1:population_size*(1-elite_ratio)
        feasible = checkFeasibility(after_intercross(z,:),s,n); 
        if (feasible == 0)          % Infeasible
             infeasible_solution = after_intercross(z,:); %disp("Infeasible Solution"); displayRouteSet(infeasible_solution,s,n); 
             feas1 = repair_infeasible_route_set(infeasible_solution, s, n, DistanceMatrix);
             if (feas1 == 0)        % Infeasibility cannot be repaired
                after_intercross(z,:) = crossover_matrix1(z,:);
             else                   % Infeasibility is repaired
                after_intercross(z,:) = feas1;
                %disp("Feasible Solution"); displayRouteSet(feas1,s,n);
             end
        end
    end

    %disp("After Inter-Crossover");
    %displayPopulation(after_intercross,n);

    end
    
    % Intra-string
    % 1. Choosing Parent population randomly
    % 2. Choosing of demarcation site randomly
    % 3. Crossover Probability P_ca 
    
    if (intracrossover == 1)
        if (intercrossover == 1) 
            crossover_matrix2 = after_intercross;
        elseif (intercrossover == 0) 
            crossover_matrix2 = current_pop_string(1:population_size,:);
        end
        
        after_intracross = zeros(population_size*(1-elite_ratio),string_length);
        
        % 1. & 3. Choosing Parent population randomly
        for w=1:population_size*(1-elite_ratio)
            r2 = rand();    %disp("r2"); disp(r2);
            if (r2 < P_ca)         % not selected as a parent
                after_intracross(w,:) = crossover_matrix2(w,:);
                %disp("Not selected as a parent"); %disp(after_intracross);
            else                   % selected as a parent
                parent_soln = crossover_matrix2(w,:); 
                %disp("Parent Solution:"); displayRouteSet(parent_soln,s,n);
                after_intracross(w,:) = intra_crossover_main(parent_soln, s, n);
                %disp("After Intra-Crossover"); displayRouteSet(after_intracross(s2,:),s,n);
            end
        end

        %disp("After Intra-Crossover"); displayPopulation(after_intracross,n);
    end
    
    % MUTATION
    pre_mutation = after_intracross;        
    after_mutation = zeros(population_size*(1-elite_ratio), s*n);
    for d=1:population_size*(1-elite_ratio)
        d1 = rand();      
        if (d1 < P_m)
            %disp("Before Mutation:"); displayRouteSet(pre_mutation(d,:),s,n);
            % mutate
            sr = randi([1,s],1);
            vec = functionRoute(pre_mutation(d,(sr-1)*n+1:sr*n));
            mp = randi([1,length(vec)],1);
            node = vec(1,mp);                               % = pre_mutation(d,(sr-1)*n+mp);
            %disp("Old Node"); disp(node);
            new_node = mutation_main(node, DistanceMatrix, vec, n);
            %disp("New Node"); disp(new_node);
            len = length(new_node);
            if (len == 1)
                after_mutation(d,:) = pre_mutation(d,:);          % other nodes stay the same
                after_mutation(d,(sr-1)*n+mp) = new_node;
            elseif (len > 1)
                after_mutation(d,:) = pre_mutation(d,:);
                after_mutation(d,(sr-1)*n+mp:(sr-1)*n+mp+len-1) = new_node;
            end
            %disp("After Mutation:"); displayRouteSet(after_mutation(d,:),s,n);
        else
            after_mutation(d,:) = pre_mutation(d,:);
        end
        %disp(after_mutation(d,:));
    end

    % checking if all nodes are in the solution set
    for z1=1:population_size*(1-elite_ratio)
        feasible = checkFeasibility(after_mutation(z1,:),s,n); 
        if (feasible == 0)                % Infeasible
             infeasible_solution = after_mutation(z1,:); 
             %disp("Infeasible Solution after Mutation"); displayRouteSet(infeasible_solution,s,n); 
             feas2 = repair_infeasible_route_set(infeasible_solution, s, n, DistanceMatrix);
             if (feas2 == 0)              % Infeasibility cannot be repaired
                after_mutation(z1,:) = pre_mutation(z1,:);
             else                         % Infeasibility is repaired
                after_mutation(z1,:) = feas2;
                %disp("Feasible Solution"); displayRouteSet(feas2,s,n);
             end
        end
    end

    %disp("After Mutation"); displayPopulation(after_mutation,n);
    
    % checking if all nodes are in the solution set
    for z1=1:population_size*(1-elite_ratio)
        j1 = ones(1,n);
        for z2=1:s*n
            for z3=1:n
                if (after_mutation(z1,z2) ~= 0)
                    if (after_mutation(z1,z2) == z3)
                        j1(1,z3) = 0;
                    end
                end
            end
        end
        %disp("j1"); disp(j1);
        if (sum(j1) ~= 0)
            after_mutation(z1,:) = pre_mutation(z1,:);
        end
    end    
    
    
    hp = population_size*elite_ratio;
    next_pop_string(hp+1:population_size,:) = after_mutation;

    %disp("New Generation"); %disp(new_generation);
    %fprintf('\n');
    %displayPopulation(next_pop_string,n);

    k = k + 1;

end


end

