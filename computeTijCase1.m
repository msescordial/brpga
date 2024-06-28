% Computing Travel Time Between Node i and Node j
% Case 1: No Transfer Needed

function [tij]=computeTijCase1(i,j,common_route, TimeMatrix) 
    
    % Initialization
    tij = Inf; 
    g = zeros(1,1);     % vector of nodes from i to j
    pos_i = 0;
    pos_j = 0;    
    k = length(common_route);
    
    % find positions of node i and j in the route
    for h=1:k
        if (common_route(1,h) == i) 
            pos_i = h; 
        end        
        if (common_route(1,h) == j) 
            pos_j = h; 
        end
    end
       
    % Case 1.1: Nodes i and j are adjacent to each other
    if ( abs(pos_i - pos_j) == 1 )
        tij = TimeMatrix(i,j);
 
    % Case 1.2: Nodes i and j are not adjacent to each other
    else
        tij = 0;
        if (pos_j > pos_i)
            g = common_route(1,pos_i:pos_j);
        else
            g = common_route(1,pos_j:pos_i);
        end        

        for m = 1: length(g)-1
            tij = tij + TimeMatrix(g(m),g(m+1));
        end
        %tij = tij + waiting_time*(length(g)-2);
        %tij = tij;
    end
    
    
end
