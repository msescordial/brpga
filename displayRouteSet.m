function displayRouteSet(solution_vector,s,n)
    
    %divide into routes
    k = 1;
    for j=1:s
        % getting the zeros
        draft_route = solution_vector(1,k:k+n-1);
        nn = nnz(draft_route); 
        C = zeros(1,nn);
        for p=1:nn
            C(1,p) = draft_route(1,p); 
        end
        for q=1:nn-1
            fprintf('%d-',C(1,q));
        end
        fprintf('%d |', C(1,nn));
        k = k + n;
    end
        fprintf('\n');
end