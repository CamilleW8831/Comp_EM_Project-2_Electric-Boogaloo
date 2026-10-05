function z = dbesselh(n, k, v)
    z = 0.5 .* (besselh(n-1, k, v) - besselh(n+1, k, v));
end