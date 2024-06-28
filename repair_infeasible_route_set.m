function [feasible_solution]=repair_infeasible_route_set(infeasible_solution, s, n, DistanceMatrix)

    feasible_solution = zeros(1,length(infeasible_solution));

    % divide the routes
    [Sc] = stringToRoutes(infeasible_solution,s,n);
    route1 = cell(s,1);
    for a=1:s
        route1{a,1}=BusRoute(Sc{a,:});
    end
    
    % determine the missing nodes
    M = zeros(1,1);        % vector of missing nodes

    g = ones(1,n);
    for b=1:n
        for c=1:s
            route2 = route1{c,:};
            d = length(route2);
            for e=1:d
                if (route2(1,e) == b)
                    g(1,b) = 0;
                end
            end
        end
    end

    %disp("g"); disp(g);

    m1 = 1;
    for f=1:n
        if (g(1,f) == 1)
            M(1,m1) = f;
            m1 = m1 + 1;
        end
    end

    %disp("Missing Nodes:"); disp(M);

    % get the neighbors of each missing node
    x = cell(1,length(M)); 

    for i=1:length(M)
        node = M(1,i);
        % vector of neighbors
        u = 1;
        vec = zeros(1,u);
        for q=1:n            
            if (DistanceMatrix(node,q) ~=0 && DistanceMatrix(node,q) ~=Inf)
                vec(1,u) = q;                  % neighbor x
                u = u + 1;
            end
        end
        x{1,i} = vec;
    end
    %disp("Neighbors:"); disp(x);

    rn = 1;
    replaced_route_no = zeros(1,rn);

    new_route_cell = cell(s,1);

    % checking each route for the neighbor x
    for i1=1:length(M)
        missing_node = M(1,i1); %disp("missing_node"); disp(missing_node);
        len = length(x{1,i1});
        for j1=1:len
            vec2 = x{1,i1};
            neighbor = vec2(j1);    %disp("neighbor"); disp(neighbor);
            % checking the routes
            for k1=1:s
                route3 = route1{k1,1};      %disp("route3"); disp(route3); 
                mem = ismember(neighbor,route3);    %disp("mem"); disp(mem); 
                if (mem == 1)
                    % Case 1 and Case 2
                    % find the position of the neighbor in that route
                    [p] = getNodePositionInARoute(route3,neighbor);   %disp("p"); disp(p);
                    lenr = length(route3);
                    new_route = zeros(1,lenr+1);
                    if (p == lenr)
                        new_route(1,1:lenr) = route3;
                        new_route(1,lenr+1) = missing_node; %disp("Case Last"); disp(new_route);
                        new_route_cell{k1,1} = new_route;
                        replaced_route_no(1,rn) = k1; rn = rn + 1;
                        break;
                    elseif (p == 1)
                        new_route(1,1) = missing_node;
                        new_route(1,2:lenr+1) = route3; %disp("Case First"); disp(new_route);
                        new_route_cell{k1,1} = new_route;
                        replaced_route_no(1,rn) = k1; rn = rn + 1;
                        break;
                    else
                    %if ((p ~= lenr) && (p ~= 1))
                        [feasible_route1] = correct_infeasible_route_case1(missing_node, neighbor, p, route3, DistanceMatrix);
                        if (feasible_route1 ~= 0)
                            new_route = feasible_route1; %disp("Case 1"); disp(new_route);
                            new_route_cell{k1,1} = new_route;
                            replaced_route_no(1,rn) = k1; rn = rn + 1;
                            break;
                        else
                            [feasible_route2] = correct_infeasible_route_case2(missing_node, neighbor, p, route3, DistanceMatrix);
                            if (feasible_route2 ~= 0)
                                new_route = feasible_route2; %disp("Case 2"); disp(new_route);
                                new_route_cell{k1,1} = new_route;
                                replaced_route_no(1,rn) = k1; rn = rn + 1;
                                break;
                            end
                        end
                    end
                end
            end
        end
    end

    replaced_route_no = sort(replaced_route_no);
    %disp("Replaced Routes"); disp(replaced_route_no);

    %disp("New Routes:"); disp(new_route_cell);

    for g2=1:s
        if (new_route_cell{g2,1} ~= 0)
            a3 = n*(g2-1)+1;            
            b3 = n*(g2-1)+length(new_route_cell{g2,1});
            feasible_solution(1,a3:b3) = new_route_cell{g2,1};
        end
    end
    
    for f2=1:s
        imr = ismember(f2,replaced_route_no);
        if (imr == 0)
            a2 = n*(f2-1)+1;            
            b2 = n*(f2-1)+n;
            feasible_solution(1,a2:b2) = infeasible_solution(1,a2:b2);
        end
    end    

    % checking if all nodes are in the feasible solution
    bin = checkFeasibility(feasible_solution,s,n);

    if (bin == 0)
        feasible_solution = 0;
    else
        %disp("Feasible Solution"); displaySolution(feasible_solution,s,n);
    end

end

