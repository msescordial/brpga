function [S1c] = stringToRoutes(S1,s,n)
    S1c = cell(s,1);
    p = 1;
    for q=1:s
        S1c{q,1} = S1(1,p:p+n-1);
        p = p + n;
    end
    %disp(S1c);
end