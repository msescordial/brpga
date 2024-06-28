function new_solution = intra_crossover_main(solution_vector, s, n)
    new_solution = zeros(1,1);
    routes_parent = zeros(s,n);         % display the routes
    p=1;
    for k=1:s
        routes_parent(k,:) = solution_vector(1,p:p+n-1);
        p = p+n;
    end
    %disp("Routes of the Parent Solution for Intra-Crossover"); 
    %disp(routes_parent);
                
    % Choose 2 routes randomly
    route_pair = 2;
    iter = 1; maxiter = 20;
             
    while (route_pair == 2)
        if (iter > maxiter)
            new_solution = solution_vector;
            break;
        end
        A = zeros(1,s); 
        for q=1:s
            A(1,q) = q;
        end
        %disp("A"); disp(A);
                
        rp = sort(randperm(numel(A),2));  %disp("rp"); disp(rp);
        x1 = A(rp(1));
        x2 = A(rp(2));

        chosen_routes = zeros(2,n);
        e = 1;
        not_chosen_routes = zeros(s-2,n);
        f = 1;
        for w=1:s
            if (w == x1)
                chosen_routes(e,:) = routes_parent(w,:);
                e = e+1;
            elseif (w == x2)
                chosen_routes(e,:) = routes_parent(w,:);
                e = e+1;
            else
                not_chosen_routes(f,:) = routes_parent(w,:);
                f = f + 1;
            end
        end
        %disp("Chosen Routes"); disp(chosen_routes);
        %disp("Not Chosen Routes"); disp(not_chosen_routes);
                
        % 2. Then Choosing of demarcation site randomly 
        % Demarcation site is at the end of the common node 
        % Only one demarcation site (for now)
        route1 = BusRoute(chosen_routes(1,:));
        route2 = BusRoute(chosen_routes(2,:));
                
        [b, h] = common_nodes(route1, route2);
        if (b == 1)
            [newroute1, newroute2] = intra_crossover_operation(route1, route2, h, n);   
            if (newroute1 == 0)
                %disp("Get another route pair.");
                route_pair = 2;
            else
                % combine the new routes into one solution route set
                C = not_chosen_routes;  %disp("C"); disp(C);
                D1 = newroute1;
                D2 = newroute2;
                route_pair = 3; 
            end
        else
            %disp("Get another route pair.");
            route_pair = 2;
        end
        iter = iter + 1;
    end

    if (iter < maxiter)
        new_solution = zeros(1,s*n);
        new_solution(1,1:n) = D1;
        new_solution(1,n+1:2*n) = D2;
        [Cr Cc] = size(C);
        f = 2;
        for q=1:Cr
            new_solution(1,f*n+1:(f+1)*n) = C(q,:); 
            f = f + 1;
        end
        %disp("New Solution After Intra-Crossover"); disp(new_solution);
        return;
    end
    if (new_solution == 0)
        new_solution = solution_vector;
    end
        
    %disp("new_solution"); disp(new_solution);
end