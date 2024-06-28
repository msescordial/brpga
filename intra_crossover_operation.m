function [newroute1, newroute2] = intra_crossover_operation(route1, route2, h, n)
    % choosing randomly 
    k = length(h);
    c = randi([1, k],1);
    for i=1:k
        if (c == i)
            demarc_node = h(i);
        end
    end
    %disp("Demarcation Node"); disp(demarc_node);
    
    m1 = length(route1); m2 = length(route2);
    % getting the position of demarc_node
    for i = 1:m1
        if (route1(i) == demarc_node)
            pos_1 = i;
        end
    end
    for j = 1:m2
        if (route2(j) == demarc_node)
            pos_2 = j;
        end
    end
      
    if ( pos_1 == m1 || pos_2 == m2)
        newroute1 = 0; newroute2 = 0;
    elseif ( pos_1 == 1 || pos_2 == 1)
        newroute1 = 0; newroute2 = 0;       
    else 
        % splitting the two routes
        part11 = route1(1,1:pos_1); p11 = length(part11);
        part12 = route1(1,pos_1+1:m1); p12 = length(part12);       
        part21 = route2(1,1:pos_2); p21 = length(part21);
        part22 = route2(1,pos_2+1:m2); p22 = length(part22);
        
        newroute1 = zeros(1,n);
        newroute2 = zeros(1,n);
        
        % intra-crossover
        newroute1(1,1:p11) = part11;
        newroute1(1,p11+1:p11+p22) = part22;
        newroute2(1,1:p21) = part21;
        newroute2(1,p21+1:p21+p12)= part12;
 
        %disp("newroute1"); disp(newroute1);
        %disp("newroute2"); disp(newroute2); 
        
        % are there two same nodes in each new route?
        A1 = BusRoute(newroute1); A2 = BusRoute(newroute2);
        A3 = sort(A1); A4 = sort(A2);      
        %disp("A3"); disp(A3); disp("A4"); disp(A4);
        
        for a5 = 1:length(A3)-1
            if (A3(1,a5) == A3(1,a5+1))
                %disp("There is a repeat node in the same route");
                newroute1 = 0; newroute2 = 0;
            end
        end
        for a6 = 1: length(A4)-1
            if (A4(1,a6) == A4(1,a6+1))
                %disp("There is a repeat node in the same route");
                newroute1 = 0; newroute2 = 0;
            end
        end  
        
    end
end