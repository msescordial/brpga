function new_node = mutation_main(node, DistanceMatrix, vec, n)
    check = zeros(1,length(vec));
    for p1=1:n
        check(1,p1) = p1;
    end
    A = ones(1,n)-ismember(check,vec);
    B = zeros(1,sum(A));
    p3 = 1;
    for p2=1:n
        if (A(1,p2) == 1)
            B(1,p3) = check(1,p2);
            p3 = p3 + 1;
        end
    end
    %disp(check); disp(A); disp(B);

    % get the nodes before and after from the vector
    for j=1:length(vec)
        if (vec(j) == node)
            pos = j;
        end
    end

    if (pos == 1)
        aft = vec(1,2);
        v1 = 1; V1 = zeros(1,v1);
        for q1=1:length(B)
            if (DistanceMatrix(B(q1),aft) ~= Inf && DistanceMatrix(B(q1),aft) ~= 0)
                V1(1,v1) = B(q1);
                v1 = v1 + 1;
            end
        end
        w1 = randi([1, length(V1)],1);
        new_node = V1(w1);
    elseif (pos == length(vec))
        bef = vec(1,length(vec)-1);
        v2 = 1; V2 = zeros(1,v2);
        for q2=1:length(B)
            if (DistanceMatrix(B(q2),bef) ~= Inf && DistanceMatrix(B(q2),bef) ~= 0)
                V2(1,v2) = B(q2);
                v2 = v2 + 1;
            end
        end
        w2 = randi([1, length(V2)],1);
        new_node = V2(w2); 
    else
        bef = vec(1,pos-1);
        aft = vec(1,pos+1);

        % checking the vicinity nodes
        f1 = 1;
        vicinity_nodes = zeros(1,f1);
        for q3=1:length(B)
            if (DistanceMatrix(bef,B(q3)) ~= Inf && DistanceMatrix(bef,B(q3)) ~= 0)
                vicinity_nodes(1,f1) = B(q3);
                f1 = f1 + 1;
            end
        end
        %disp("Vicinity Nodes"); disp(vicinity_nodes);

        if (vicinity_nodes == 0)
            new_node = node;
            return;
        end

        len = length(vicinity_nodes);
    
        f2 = 1;
        candidates = zeros(1,f2);

        % Check if a single link exists
        for k=1:len
            if (DistanceMatrix(vicinity_nodes(1,k),aft) ~= Inf && DistanceMatrix(vicinity_nodes(1,k),aft) ~= 0)
                candidates(1,f2) = vicinity_nodes(1,k);
                f2 = f2 + 1;
            end
        end
        if (candidates ~= 0)
            r = randi([1, length(candidates)],1);
            new_node = candidates(r);
            return;
        end
        %disp("Candidates:"); disp(candidates);

        % If no single link,
        b = randi([1, len], 1);
        %disp("b"); disp(b);
        vici = vicinity_nodes(b);    
        f3 = 1;      
        k_paths = 5;
        [shortestPaths, totalCosts] = kSP(DistanceMatrix, vici, aft, k_paths);
        for h=1:k_paths
            path = shortestPaths{1,k}; len_path = length(path);
            %disp("path"); disp(path);
            n_node = path(1,1:len_path-1); 
            D = ismember(n_node,vec);
            if (sum(D) == 0)
                new_node = n_node;
                return;
            else 
                new_node = node;
                return;
            end
        end
    end
    
    if (new_node == 0)
        new_node = node;
    end

    %disp("New Node"); disp(new_node);
end

