% Note: route must have nonzero elements 

function [feasible_route1] = correct_infeasible_route_case1(missing_node, neighbor, p, route, DistanceMatrix)
    len = length(route);
    feasible_route1 = zeros(1,len+1);

    %disp("route"); disp(route);

    s1 = route(1,1:p-1);
    s2 = route(1,p+1:len);

    reversed_s1 = zeros(1,p-1);
    appended_s2 = zeros(1,len-p);

    for i=1:p-1
        reversed_s1(1,i) = route(1,p-i);
    end
    %disp("S1 and its reverse"); disp(s1); disp(reversed_s1);
    
    neighbor_node = route(1,1);
    for j=p+1:len
        ap_node = route(1,j);
        if ( DistanceMatrix(ap_node, neighbor_node) ~= 0 && DistanceMatrix(ap_node, neighbor_node) ~= Inf)
            appended_s2(1,j-p) = ap_node;
            ap_node = neighbor_node;
        end
    end
    %disp("S2 and appended S2"); disp(s2); disp(appended_s2);

    if (sum(appended_s2) == 0)
        feasible_route1 = 0;
    else 
        % Assignment
        feasible_route1(1,1) = missing_node;
        feasible_route1(1,2) = neighbor;
        feasible_route1(1,3:p+1) = reversed_s1;
        feasible_route1(1,p+2:len+1) = appended_s2;
    end

    %disp("feasible_route1"); disp(feasible_route1); 

end