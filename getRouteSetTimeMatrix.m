function SolutionTimeMatrix = getRouteSetTimeMatrix(S0r,s,TimeMatrix, transfer_time)
    n = size(TimeMatrix,1);
        
    % Initialization
    SolutionTimeMatrix = Inf([n,n]);
    %disp("Initialization"); disp(SolutionTimeMatrix);
    
    % Listing the routes in S0 in Matrix Form
    B = zeros(s,n);
    for t=1:s
        B(t,:) = S0r{t,1};        % all s bus routes
    end
    % fprintf('B: \n'); disp(B);
    
    
    % Determining the Shortest Time from Node i to Node j

    % 1. Diagonals must be zero
    for i=1:n
    for j=1:n
        if (i == j)
        	SolutionTimeMatrix(i,j) = 0; 
        end
    end
    end

    %disp("SolutionTimeMatrix After Case 1");
    %disp(SolutionTimeMatrix); 
    
    % 2. Identifying Routes which contain Node i and Node j, respectively 
    %    and Case 1: No Transfer is Needed
    for i=1:n
        for j=1:n           
            % Identifying the Routes where Nodes i and j are in  
            f1 = 1;
            f2 = 1;
            rin = zeros(1,f1);       % vector which stores the routei nos.
            rjn = zeros(1,f2);       % vector which stores the routej nos.
            for p=1:s
                routep = functionRoute(B(p,:));  
                mi = ismember(routep,i);
                mj = ismember(routep,j);                   
                if (sum(mi) == 1)       % node i is in route p
                    rin(1,f1) = p;
                    f1 = f1+1;
                end
                if (sum(mj) == 1)       % node j is in route p
                    rjn(1,f2) = p;
                    f2 = f2+1;
                end
            end   
            %disp("i"); disp(i); disp("j"); disp(j); 
            %disp("Routes containing i"); disp(rin); disp("Routes containing j"); disp(rjn);
            
            % Case 1: No Transfer is Needed
            if (sum(ismember(rin,rjn)) > 0)     % common route         
            	f3 = 1;
            	crn = zeros(2,f3);
            	for u = 1:length(rin)
                for v = 1: length(rjn)
                    if (rin(u) == rjn(v))
                    	crn(1,f3) = rin(u); 
                    	f3 = f3 + 1;
                    end
                end
                end      
                [r1 c1] = size(crn);
                for w = 1:c1
                    p1 = crn(1,w);
                    common_route = functionRoute(B(p1,:));
                    %crn(2,w) = computeTijCase1(i,j,common_route, waiting_time, TimeMatrix);
                    crn(2,w) = computeTijCase1(i,j,common_route, TimeMatrix);
                end
                %disp("Case 1: No Transfer is Needed");
                %disp("crn"); disp(crn);
                SolutionTimeMatrix(i,j) = min(crn(2,:));
            
            elseif (sum(ismember(rin,rjn)) == 0)            % there is no common route                        
            	f4 = 1;
            	crn2 = zeros(3,f4);
                for u = 1: length(rin)
                for v = 1: length(rjn)
                	pi = rin(1,u); 
                    pj = rjn(1,v); 
                    routei = functionRoute(B(pi,:)); 
                    routej = functionRoute(B(pj,:));                            
                    %disp("Are there common nodes?");disp(sum(ismember(routei,routej)));
                    if (sum(ismember(routei,routej)) > 0)   % there is at least one common node
                        % Case 2: One Transfer is Needed
                        %disp("Case 2: One Transfer is Needed"); 
                        crn2(1,f4) = pi;
                        crn2(2,f4) = pj;
                        %crn2(3,f4) = computeTijCase2(i,j, routei, routej, waiting_time, transfer_time, TimeMatrix);
                        crn2(3,f4) = computeTijCase2(i,j, routei, routej, transfer_time, TimeMatrix);
                        f4 = f4+1;
                        %disp("crn2"); disp(crn2);
                        SolutionTimeMatrix(i,j) = min(crn2(3,:)); 
                    end       
                end
                end                 
            end
    	%end
        end
    end
    
    %disp("SolutionTimeMatrix After Cases 1 and 2");
    %disp(SolutionTimeMatrix); 
    
        
    % 3. Case 3: Two Transfers Are Needed
    for i=1:n
    for j=1:n
         if ( abs(i-j) > 0 ) 
         if (SolutionTimeMatrix(i,j) == Inf)
            %disp("i"); disp(i); disp("j"); disp(j);
            %disp("Case 3?");
            f5 = 1;
            f6 = 1;
            rin = zeros(1,f5);       % vector which stores the routei nos.
            rjn = zeros(1,f6);       % vector which stores the routej nos.
            for p1=1:s
                routep = functionRoute(B(p1,:));  
                mi = ismember(routep,i);
                mj = ismember(routep,j);                   
                if (sum(mi) == 1)       % node i is in route p
                    rin(1,f5) = p1;
                    f5 = f5+1;
                end
                if (sum(mj) == 1)       % node j is in route p
                    rjn(1,f6) = p1;
                    f6 = f6+1;
                end
            end   
            %disp("i"); disp(i); disp("j"); disp(j); 
            %disp("Routes containing i"); disp(rin); disp("Routes containing j"); disp(rjn);
            
            for u3 = 1: length(rin)
            for v3 = 1: length(rjn)
            	pi3 = rin(1,u3); 
                pj3 = rjn(1,v3);
                routei = functionRoute(B(pi3,:));
                routej = functionRoute(B(pj3,:)); 
                
                %disp("Case 3: Two Transfers are Needed");
                
                % getting the possible transfer routes
                % (excluding pi3 and pj3)
                f8 = 1;
                tr = zeros(1,f8);
                for g=1:s
                    if ( abs(g - pi3) > 0 && abs(g - pj3) > 0)
                        tr(1,f8) = g;
                        f8 = f8 + 1;
                    end
                end              
                %disp("Possible Transfer Routes"); disp(tr);
                
                f7=1;
                crn3 = zeros(4,f7);
                
                for g1 = 1: length(tr)
                	ptf = tr(g1);
                    routetf = functionRoute(B(ptf,:));
                    if (sum(ismember(routei,routetf))>0 && sum(ismember(routej,routetf))>0) 
                        routef = routetf;
                        crn3(1,f7) = pi3;
                        crn3(2,f7) = ptf;
                        crn3(3,f7) = pj3;
                        %crn3(4,f7) = computeTijCase3(i,j,routei, routef, routej, waiting_time, transfer_time, TimeMatrix);
                        crn3(4,f7) = computeTijCase3(i,j,routei, routef, routej, transfer_time, TimeMatrix);
                        f7 = f7 + 1;
                    end                               
                end
                %disp("crn3"); disp(crn3);
                SolutionTimeMatrix(i,j) = min(crn3(4,:));
            end
            end
                   
         end            
         end
    end
    end
    
    for i=1:n
    for j=1:n
         if ( abs(i-j) > 0 ) 
            if (SolutionTimeMatrix(i,j) == 0)
                SolutionTimeMatrix(i,j) = Inf;
            end
         end
    end
    end

    %disp("SolutionTimeMatrix After Case 3");
    %disp(SolutionTimeMatrix);      
end