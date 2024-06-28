function [tij]=computeTijCase3(i,j,routei, routef, routej, transfer_time, TimeMatrix)

    % Initialization
    tij = Inf; 
    
    v1 = 1;
    h1 = zeros(1,1);     % vector of common nodes of routei and routef
    
    v2 = 1;
    h2 = zeros(1,1);     % vector of common nodes of routef and routej
    
    for p1 = 1:length(routei)
        for q1 = 1:length(routef)
            if (routei(p1) == routef(q1))
                h1(1,v1) = routei(1,p1);     % common node 
                v1 = v1+1;
            end
        end
    end
    
    for p2 = 1:length(routef)
        for q2 = 1:length(routej)
            if (routef(p2) == routej(q2))
                h2(1,v2) = routef(1,p2);     % common node 
                v2 = v2+1;
            end
        end
    end
   
    k1 = length(h1);  
    k2 = length(h2); 
            
    % Matrix M of travel times
    % 1st and 2nd row: node h1 and h2   (all possible combinations)
    % 3rd row: travel time between nodes i and h1
    % 4th row: travel time between nodes h1 and h2
    % 5th row: travel time between nodes h2 and j
    % 6th row: Total travel time between nodes i and j
    
    M = zeros(6,k1*k2);  

    
    % 1st and 2nd row of M:
    w = 1;
    for u = 1:k1
        for v = w:k2+w-1
            M(1,v)=h1(u);
        end
        for y = 1:k2       
            M(2,y+w-1)=h2(y);
        end
        w = w + k2;
    end
      
    % 3rd row of M  
    for a=1:k1*k2
        [tih1]=tijCase1(i,M(1,a),routei, TimeMatrix);
        M(3,a)=tih1;
    end 
    
    % 4th row of M
    for b=1:k1*k2
        [th1h2]=tijCase1(M(1,b),M(2,b),routef, TimeMatrix);
        M(4,b)=th1h2;
    end
    
    % 5th row of M
    for c=1:k1*k2
        [th2j]=tijCase1(M(2,c),j,routej, TimeMatrix);
        M(5,c)=th2j;
    end
    
    % 6th row of M
    for d=1:k1*k2
        M(6,d)=M(3,d)+ transfer_time + M(4,d) + transfer_time + M(5,d);
    end   
    
    %disp("i"); disp("j"); disp(i); disp(j);
    %disp("M"); disp(M); 
    tij = min(M(6,:)); %disp("tij"); disp(tij);
end