function displayInitialRouteSet(S0,BusRouteID)
    s = size(S0,1);
    for b = 1:s
        r = S0(b,1);
        draft_route = BusRouteID{r,2};
        
        fprintf('Route %d:', r);
        displayRoute(draft_route);
    end
end