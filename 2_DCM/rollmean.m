function res = rollmean(v, k)

    v = v(~isnan(v));
    n = length(v);
    
    if n < 1, res = 0; return; end
    if n < k, res = NaN; return; end
    if n == 1, res = v; return; end
    
    res = zeros(n-k+1, 1);
    for i = k:n
        res(i-k+1) = mean( v(i-k+1:i) );
    end         

end
