% Note: route must have nonzero elements 

function [feasible_route2] = correct_infeasible_route_case2(missing_node, neighbor, p, route, DistanceMatrix)
    len = length(route);
    feasible_route2 = zeros(1,len+1);

    %disp("route"); disp(route);

    s1 = route(1,1:p-1);
    s2 = route(1,p+1:len);

    reversed_s2 = zeros(1,len-p);
    appended_s1 = zeros(1,p-1);

    for i=1:len-p
        reversed_s2(1,i) = route(1,len+1-i);
    end
    %disp("S2 and its reverse"); disp(s2); disp(reversed_s2);
    
    neighbor_node = route(1,len);
    for j=1:p-1
        ap_node = route(1,j);
        if ( DistanceMatrix(ap_node, neighbor_node) ~= 0 && DistanceMatrix(ap_node, neighbor_node) ~= Inf)
            appended_s1(1,p-j) = ap_node;
            ap_node = neighbor_node;
        end
    end
    %disp("S1 and appended S1"); disp(s1); disp(appended_s1);

    if (sum(appended_s1) == 0)
        feasible_route2 = 0;
    else 
        % Assignment
        feasible_route2(1,1:p-1) = appended_s1;
        feasible_route2(1,p:len-1) = reversed_s2;
        feasible_route2(1,len) = neighbor;
        feasible_route2(1,len+1) = missing_node;
    end

    %disp("feasible_route1"); disp(feasible_route2); 

end