function [DistanceMatrix,TimeMatrix,TravelDemandMatrix,TerminalNodes,k,s,transfer_time]=network_manila_2()

%------------------------------ Remarks: ----------------------------------
% Network Inputs: 
%  d = Stop-to-Stop Distance Matrix
%  T = Stop-to-Stop Time Matrix
%  D = Travel Demand Matrix
%  t = Terminal Nodes
%  k = for k-shortest paths algorithm: this refers to the first k
%    shortest paths
%  s = no. of bus routes for the network
%  transfer_time = time for each transfer (constant)

% T
load("final_time.mat","-mat");
TimeMatrix = array;
for i=1:59
    for j=1:59
        if abs(i-j)>0
            if TimeMatrix(i,j) == 0
                TimeMatrix(i,j) = inf;
            end
        end
    end
end
disp("Time Matrix"); disp(TimeMatrix);


% d
load("final_distance.mat","-mat");
DistanceMatrix = array;
for i=1:59
    for j=1:59
        if abs(i-j)>0
            if DistanceMatrix(i,j) == 0
                DistanceMatrix(i,j) = inf;
            end
        end
    end
end
disp("Distance Matrix"); disp(DistanceMatrix);

% D 
load("final_demand.mat","-mat");
TravelDemandMatrix = array;
disp("Travel Demand Matrix"); disp(TravelDemandMatrix);

% t
TerminalNodes = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];      % terminal nodes
%TerminalNodes = [2 3 4 5 7 8 9 11 12 13 14 16 17 18 19 20 21 24 25];

k = 4;                  % k shortest Paths for each node to node
s = 20;                 % no. of routes in a bus network
transfer_time = 5;      % transfer time is 5 minutes

end


