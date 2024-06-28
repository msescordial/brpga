function [S1] = GeneticAlgo(population_size, S0, S0r, maxiter, n, DistanceMatrix, TravelDemandMatrix, TimeMatrix, BusRouteID, TotalNoOfRoutes, s, transfer_time)
    s = length(S0); 
    
    init_pop_matrix = cell(population_size,3);
    
    % Creating the Initial Population 
    % 1st col: no., 2nd col: pop, 3rd col: Obj Func Value, 4th col: Fitness Value
    % First Row
    init_pop_matrix{1,1} = 1;
    init_pop_matrix{1,2} = transpose(S0);
    S0_SolutionTimeMatrix = TotalTime(S0r,s,TimeMatrix, transfer_time);
    ofv0 = ObjFuncVal(S0r,TravelDemandMatrix,TimeMatrix,S0_SolutionTimeMatrix,n);
    init_pop_matrix{1,3} = ofv0;   init_pop_matrix{1,4} = ofv0;   
    % Succeeding Rows
    for g = 1: population_size-1
        init_pop_matrix{g+1,1} = g+1;
        [S1,S1r] = generateInitialNetwork(DistanceMatrix, BusRouteID, TotalNoOfRoutes,s);
        init_pop_matrix{g+1,2} = transpose(S1);
        S1_SolutionTimeMatrix = TotalTime(S1r,s,TimeMatrix, transfer_time);
        ofv1 = ObjFuncVal(S1r,TravelDemandMatrix,TimeMatrix,S1_SolutionTimeMatrix,n);
        init_pop_matrix{g+1,3} = ofv1; init_pop_matrix{g+1,4} = ofv1;  
    end
    %disp("Initial Population Matrix"); disp(init_pop_matrix);
    
    % SELECTION
    % Rank1 = Highest Ftness Value
    mx = zeros(1,population_size);
    my = cell(1,population_size);
    for i=1:population_size
        mx(1,i) = init_pop_matrix{i,4};
        my{1,i} = init_pop_matrix{i,2};
    end
    [mxsorted,I] = sort(mx);
    mysorted = my(I);
    %disp("mxsorted"); disp(mxsorted); %disp("mysorted"); disp(mysorted);       
    parent_pop = zeros(population_size*0.5, s);
    for j=1:population_size*0.5
        parent_pop(j,:)=mysorted{1,j};
    end
    %disp("Parent Population:"); disp(parent_pop);
    
    % SOLUTION REPRESENTATION FOR GA 
    len_string = s*n;
    solution_s0 = zeros(population_size,len_string);
    for h=1:population_size*0.5
        y1=1; y2=n;
        for k=1:s
            route_string = BusRouteID{parent_pop(h,k),2};
            solution_s0(h,y1:y2)=route_string;
            y1=y1+n; y2=y2+n;
        end
    end
    %disp("Parent Population String"); disp(solution_s0);
    % Note: Demarcation Line at the end of every nth node

    iter = 1;

    % GA LOOP
    while (iter <= maxiter)
    disp("Iteration:"); disp(iter);
    if (iter ~= 1)
        % compute  for the objective function value, and sort
        gen_routes = cell(population_size,2);
        for ov3=1:population_size
            gen_routes{ov3,1} = new_generation(ov3,:);
        end
        for ov=1:population_size
            sor = cell(s,1);
            p_s = 1;
            for ov2=1:s
                sor{ov2,1} = new_generation(ov,p_s:p_s+n-1);
                p_s = p_s+n;
            end
            %disp("route sor"); disp(sor);
            sor_SolutionTimeMatrix = TotalTime(sor,s,TimeMatrix, transfer_time);
            %disp("Solution Time Matrix"); disp(sor_SolutionTimeMatrix);
            obj_val = ObjFuncVal(sor,TravelDemandMatrix,DistanceMatrix,sor_SolutionTimeMatrix,n);
            gen_routes{ov,2} = obj_val;
        end
        % Rank1 = Highest Fitness Value
        rx = zeros(1,population_size);
        ry = cell(1,population_size);
        for i=1:population_size
            rx(1,i) = gen_routes{i,2};
            ry{1,i} = gen_routes{i,1};
        end
        [rxsorted,I] = sort(rx);
        rysorted = ry(I);
        %disp("rxsorted"); disp(rxsorted); disp("rysorted"); disp(rysorted);       
        
        sorted_pop = zeros(population_size*0.5, s*n);
        for j=1:population_size*0.5
            sorted_pop(j,:)=rysorted{1,j};
        end
        %disp("Sorted Population:"); disp(sorted_pop);
        solution_s0 = sorted_pop;
        intercrossover = 1;

        % display the best solution for each iteration/generation
        disp("Best Solution for this Iteration"); %displaySolution(sorted_pop(1,:),s,n);
        disp(rxsorted(1,1));

        %diff = var(rxsorted); disp(diff);
        %if (diff < 10^(-10))
        %    break;
        %end

    elseif (iter == 1) 
        % ready for crossover 
        intercrossover = 1;
    end
     

    % CROSSOVER
    % Inter-string
    % 1. Rearrangement of Parent population randomly
    % 2. Choosing of demarcation site randomly
    % 3. Crossover Probability P_ce 
    
    if (intercrossover == 1)
    % 1.   
    crossover_matrix = zeros(population_size*0.5,len_string);
    r = randperm(population_size*0.5);      
    for m=1:population_size*0.5
        crossover_matrix(m,:) = solution_s0(r(m),:);
    end
    %disp("Crossover Rearrangement"); %disp(crossover_matrix);
    %displayPopulationSolution(crossover_matrix,n);
    
    % 2. & 3.
    % Demarcation site is at the end of nth element, 2n, 3n, m*(s-1)
    % Only one demarcation site (for now)
    P_ce = 0.5;     % predefined
    
    % Crossover for each pair
    after_intercross = zeros(population_size*0.5,len_string);
    for s1=1:population_size*0.5*0.5
        crossover_strings = zeros(2,s*n);
        crossover_strings(1,:)=crossover_matrix(2*s1-1,:);
        crossover_strings(2,:)=crossover_matrix(2*s1,:);
        %disp("Crossover Strings"); displaySolution(crossover_strings(1,:),s,n); displaySolution(crossover_strings(2,:),s,n);
        %r = rand; disp("r"); disp(r);
        if ( r > P_ce )         % then crossover 
            new_string = zeros(2,s*n);
            demarc_line = n*randi([1,s-1],1) + 1;
            %disp("Demarcation Line"); disp(demarc_line);
            new_string(1,1:demarc_line-1) = crossover_strings(1,1:demarc_line-1);
            new_string(1,demarc_line:s*n) = crossover_strings(2,demarc_line: s*n);
            new_string(2,1:demarc_line-1) = crossover_strings(2,1:demarc_line-1);
            new_string(2,demarc_line:s*n) = crossover_strings(1,demarc_line: s*n);
                    
            %disp("Crossover String"); disp(crossover_strings);
            %disp("New String"); disp(new_string);
            
            after_intercross(2*s1-1,:) = new_string(1,:);
            after_intercross(2*s1,:) = new_string(2,:);            
        else
            after_intercross(2*s1-1,:) = crossover_strings(1,:);
            after_intercross(2*s1,:) = crossover_strings(2,:);
        end
        %disp("After Inter-Crossover");
        %displaySolution(after_intercross(2*s1-1,:),s,n);
        %displaySolution(after_intercross(2*s1,:),s,n);
    end

    % checking if all nodes are in the solution set
    for z=1:population_size*0.5
        bin = checkAllNodes(after_intercross(z,:),s,n); %disp("bin"); disp(bin);
        if (bin == 0)
             infeasible_solution = after_intercross(z,:); 
             %disp("Infeasible Solution"); displaySolution(infeasible_solution,s,n); 
             Aic = repair_infeasibility(infeasible_solution, s, n, DistanceMatrix);
             if (Aic == 0)
                after_intercross(z,:) = crossover_matrix(z,:);
             else
                after_intercross(z,:) = Aic;
                %disp("Feasible Solution"); displaySolution(Aic,s,n);
             end
        end
    end

    %disp("After Inter-Crossover"); displayPopulationSolution(after_intercross,n);

    end
    
    % Intra-string
    % 1. Choosing Parent population randomly
    % 2. Choosing of demarcation site randomly
    % 3. Crossover Probability P_ca 
    
    intracrossover = 1;
    P_ca = 0.5;
    
    if (intracrossover == 1)
        if (intercrossover == 1) 
            crossover_matrix2 = after_intercross;
        else
            crossover_matrix2 = solution_s0(1:population_size*0.5,:);
        end
        
        after_intracross = zeros(population_size*0.5,len_string);
        
        % 1. & 3. Choosing Parent population randomly
        for s2=1:population_size*0.5
            r2 = rand(); %disp("r2"); disp(r2);
            if (r2 < P_ca) 
                % not selected as a parent
                after_intracross(s2,:) = crossover_matrix2(s2,:);
                %disp("Not selected as a parent"); %disp(after_intracross);
            else 
                % selected as a parent
                parent_soln = crossover_matrix2(s2,:); 
                %disp("Parent Solution:"); displaySolution(parent_soln,s,n);
                %disp("r"); disp(r2);
                after_intracross(s2,:) = intra_crossover_main(parent_soln, s, n);
                %disp("After Intra-Crossover"); displaySolution(after_intracross(s2,:),s,n);
            end
        end

        %disp("After Intra-Crossover"); displayPopulationSolution(after_intracross,n);
    end
    
    % MUTATION
    pre_mutation = after_intracross;        
    after_mutation = zeros(population_size*0.5, s*n);
    P_m = 0.05;            % mutation probability
    for d=1:population_size*0.5
        d1 = rand();      
        if (d1 < P_m)
            %disp("Before Mutation:"); displaySolution(pre_mutation(d,:),s,n);
            % mutate
            sr = randi([1,s],1);
            vec = BusRoute(pre_mutation(d,(sr-1)*n+1:sr*n));
            mp = randi([1,length(vec)],1);
            node = vec(1,mp);                               % = pre_mutation(d,(sr-1)*n+mp);
            %disp("Old Node"); disp(node);
            new_node = mutation_op(node, DistanceMatrix, vec, n);
            %disp("New Node"); disp(new_node);
            len = length(new_node);
            if (len == 1)
                after_mutation(d,:) = pre_mutation(d,:);          % other nodes stay the same
                after_mutation(d,(sr-1)*n+mp) = new_node;
            elseif (len > 1)
                after_mutation(d,:) = pre_mutation(d,:);
                after_mutation(d,(sr-1)*n+mp:(sr-1)*n+mp+len-1) = new_node;
            end
            %disp("After Mutation:"); displaySolution(after_mutation(d,:),s,n);
        else
            after_mutation(d,:) = pre_mutation(d,:);
        end
        %disp(after_mutation(d,:));
    end

    % checking if all nodes are in the solution set
    for z1=1:population_size*0.5
        bin = checkAllNodes(after_mutation(z1,:),s,n); %disp("bin"); disp(bin);
        if (bin == 0)
             infeasible_solution = after_mutation(z1,:); 
             %disp("Infeasible Solution after Mutation"); displaySolution(infeasible_solution,s,n); 
             Aic2 = repair_infeasibility(infeasible_solution, s, n, DistanceMatrix);
             if (Aic2 == 0)
                after_mutation(z1,:) = pre_mutation(z1,:);
             else
                after_mutation(z1,:) = Aic2;
                %disp("Feasible Solution"); displaySolution(Aic2,s,n);
             end
        end
    end

    %disp("After Mutation"); displayPopulationSolution(after_mutation,n);

    new_generation = solution_s0;
    hp = population_size*0.5;
    
    % checking if all nodes are in the solution set
    for z1=1:population_size*0.5
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

    new_generation(hp+1:hp*2,:) = after_mutation;
    %disp("New Generation"); %disp(new_generation);
    displayPopulationSolution(new_generation,n);


    iter = iter + 1;

    end

    % Choose the best solution in the new_generation
        % compute  for the objective function value, and get the maximum
        gen_routes = cell(population_size,2);
        for ov3=1:population_size
            gen_routes{ov3,1} = new_generation(ov3,:);
        end
        for ov=1:population_size
            sor = cell(s,1);
            p_s = 1;
            for ov2=1:s
                sor{ov2,1} = new_generation(ov,p_s:p_s+n-1);
                p_s = p_s+n;
            end
            %disp("route sor"); disp(sor);
            sor_SolutionTimeMatrix = TotalTime(sor,s,TimeMatrix, transfer_time);
            %disp("Solution Time Matrix"); disp(sor_SolutionTimeMatrix);
            obj_val = ObjFuncVal(sor,TravelDemandMatrix,DistanceMatrix,sor_SolutionTimeMatrix,n);
            gen_routes{ov,2} = obj_val;
        end
        % Rank1 = Highest Fitness Value
        rx = zeros(1,population_size);
        ry = cell(1,population_size);
        for i=1:population_size
            rx(1,i) = gen_routes{i,2};
            ry{1,i} = gen_routes{i,1};
        end
        [rxsorted,I] = sort(rx);
        rysorted = ry(I);
        %disp("rxsorted"); disp(rxsorted); disp("rysorted"); disp(rysorted);       
        
        sorted_pop = zeros(population_size, s*n);
        for j=1:population_size
            sorted_pop(j,:)=rysorted{1,j};
        end
        %disp("Sorted Population:"); %disp(sorted_pop);
        %displayPopulationSolution(sorted_pop,n);

    S1 = sorted_pop(1,:);
    %disp("S1"); disp(S1);

end

