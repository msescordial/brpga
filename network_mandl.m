function [DistanceMatrix,TimeMatrix,TravelDemandMatrix,TerminalNodes,k,s,transfer_time]=network_mandl()

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
TimeMatrix = [0 8 inf inf inf inf inf inf inf inf inf inf inf inf inf;
                    8 0 2 3 6 inf inf inf inf inf inf inf inf inf inf;
                    inf 2 0 inf inf 3 inf inf inf inf inf inf inf inf inf;
                    inf 3 inf 0 4 4 inf inf inf inf inf 10 inf inf inf;
                    inf 6 inf 4 0 inf inf inf inf inf inf inf inf inf inf;
                    inf inf 3 4 inf 0 inf 9 inf inf inf inf inf inf 3;
                    inf inf inf inf inf inf 0 inf inf 7 inf inf inf inf 2;
                    inf inf inf inf inf 9 inf 0 inf 8 inf inf inf inf 2;
                    inf inf inf inf inf inf inf inf 0 inf inf inf inf inf 8;
                    inf inf inf inf inf inf 7 8 inf 0 5 inf 10 8 inf;
                    inf inf inf inf inf inf inf inf inf 5 0 10 5 inf inf;
                    inf inf inf 10 inf inf inf inf inf inf 10 0 inf inf inf;
                    inf inf inf inf inf inf inf inf inf 10 5 inf 0 2 inf;
                    inf inf inf inf inf inf inf inf inf 8 inf inf 2 0 inf;
                    inf inf inf inf inf 3 2 2 8 inf inf inf inf inf 0];
% d
DistanceMatrix = 0.3*TimeMatrix;        % base on the average travel speed of 18 kph

% D
TD = [0 400 200 60 80 150 75 75 30 160 30 25 35 0 0;
                      400 0 50 120 20 180 90 90 15 130 20 10 10 5 0;
                      200 50 0 40 60 180 90 90 15 45 20 10 10 5 0;
                      60 120 40 0 50 100 50 50 15 240 40 25 10 5 0;
                      80 20 60 50 0 50 25 25 10 120 20 15 5 0 0;
                      150 180 180 100 50 0 100 100 30 880 60 15 15 10 0;
                      75 90 90 50 25 100 0 50 15 440 35 10 10 5 0;
                      75 90 90 50 25 100 50 0 15 440 35 10 10 5 0;
                      30 15 15 15 10 30 15 15 0 140 20 5 0 0 0;
                      160 130 45 240 120 880 440 440 140 0 600 250 500 200 0;
                      30 20 20 40 20 60 35 35 20 600 0 75 95 15 0;
                      25 10 10 25 15 15 10 10 5 250 75 0 70 0 0;
                      35 10 10 10 5 15 10 10 0 500 95 70 0 45 0;
                      0 5 5 5 0 10 5 5 0 200 15 0 45 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

TravelDemandMatrix = TD*0.20;

% t
TerminalNodes = [1 8 9 12 14];      % terminal nodes


k = 4;                  % k shortest Paths for each node to node
s = 4;                  % no. of routes in a bus network
transfer_time = 5;      % transfer time is 5 minutes

end


