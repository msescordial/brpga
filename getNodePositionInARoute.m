function [p] = getNodePositionInARoute(route,node)
    len = length(route);
    
    for i=1:len
        if (route(1,i) == node)
            p = i;
        end
    end

end