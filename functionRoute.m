function [final_route]=functionRoute(draft_route)
    % Input is a row vector, m must be 1
    [m n] = size(draft_route);
    %fprintf('Size is %d by %d: \n',m, n);
    if (m == 1)
        nn = nnz(draft_route); 
        final_route = zeros(1,nn);
        for i=1:nn
            final_route(1,i) = draft_route(1,i); 
        end
    end
end
