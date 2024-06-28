function displayRoute(draft_route)
    % Input is a row vector, m must be 1
    [m n] = size(draft_route);

    if (m == 1)
        nn = nnz(draft_route); 
        C = zeros(1,nn);
        for i=1:nn
            C(1,i) = draft_route(1,i); 
        end
        % disp(C);
        for i=1:nn-1
            fprintf('%d-',C(1,i));
        end
        fprintf('%d\n',C(1,nn));
    end
end
