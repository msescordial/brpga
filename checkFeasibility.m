function bin = checkFeasibility(solution_vector,s,n)

        % checking if all nodes are in the solution vector
        nodes = ones(1,n);
        for a=1:s*n
        for b=1:n
            if (solution_vector(1,a) ~= 0)
            if (solution_vector(1,a) == b)
                nodes(1,b) = 0;
            end
            end
        end
        end

        if (sum(nodes) ~= 0)
            bin = 0;
        else
            bin = 1;
        end


end