% Computing Travel Time Between Node i and Node j
% Case 2: One Transfer is Needed

function [tij]=computeTijCase2(i,j, routei, routej, transfer_time, TimeMatrix) 
    
    % Initialization
    tij = Inf; 
    
    v = 1;
    h = zeros(1,1);     % vector of common nodes of routei and routej
    
    for p = 1:length(routei)
        for q = 1:length(routej)
            if (routei(p) == routej(q))
                h(1,v) = routei(1,p);     % common node 
                v = v+1;
            end
        end
    end
    
    k = length(h);    
    %disp("h"); disp(h);
            
    % Matrix M of travel times
    % 1st row: travel time between nodes i and h
    % 2nd row: travel time between nodes h and j
    % 3rd row: Total travel time between nodes i and j
    M = zeros(3,k);     
    
    % 1st row of M  
    for a=1:length(h)
        %[tih]=tijCase1(i,h(a),routei, waiting_time, TimeMatrix);
        [tih]=tijCase1(i,h(a),routei, TimeMatrix);
        M(1,a)=tih;
    end
    
    % 2nd row of M
    for b=1:length(h)
        %[thj]=tijCase1(h(b),j,routej, waiting_time, TimeMatrix);
        [thj]=tijCase1(h(b),j,routej, TimeMatrix);
        M(2,b)=thj;
    end
    
    % 3rd row of M
    for c=1:length(h)
        M(3,c)=M(1,c)+ transfer_time + M(2,c);
    end

    %disp("i"); disp("j"); disp(i); disp(j);
    %disp("M"); disp(M);    
    tij = min(M(3,:)); %disp("tij"); disp(tij);
end
