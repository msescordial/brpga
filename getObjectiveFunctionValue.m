function E0 = getObjectiveFunctionValue(S0r,TravelDemandMatrix,DistanceMatrix,SolutionTimeMatrix,n)
    cp = 1;       % unit cost of passenger travel time per minute (Php/pass.-min)
    cb = 5;       % unit bus operating cost per km (Php/bus-km)
    
    % Total Passenger Cost
    E1 = 0;
    for i=1:n
        for j=1:n
            E1 = E1 + TravelDemandMatrix(i,j)*(SolutionTimeMatrix(i,j)/60);
            %E1 = E1 + TravelDemandMatrix(i,j)*(SolutionTimeMatrix(i,j));
        end
    end
    E1 = cp*E1;
    
    % Total Operating Cost
    E2 = 0;
    % Listing the routes in S0 in Matrix Form
    [s c] = size(S0r);
    B = zeros(s,n);
    for t=1:s
        B(t,:) = S0r{t,1};        % all s bus routes
    end
    %fprintf('B: \n'); disp(B);
    
    dist = zeros(s,1);
    for t2=1:s
        fin_route = functionRoute(B(t2,:));
        for t3=1:length(fin_route)-1
            dist(t2,1) = dist(t2,1) + DistanceMatrix(fin_route(t3),fin_route(t3+1));
        end
    end
    %fprintf('Distance Travelled by Routes: \n'); disp(dist);    
    E2 = cb*(2*sum(dist))/1000;
    %E2 = cb*(2*sum(dist));
    
    E0 = E1 + E2;



    % Average Travel Time 
%    E1 = 0;
%    for i=1:n
%        for j=1:n
%            E1 = E1 + TravelDemandMatrix(i,j)*SolutionTimeMatrix(i,j);
%        end
%    end
%    TotalDemand = sum(TravelDemandMatrix, "all");
%    E0 = E1/TotalDemand;
    
end
