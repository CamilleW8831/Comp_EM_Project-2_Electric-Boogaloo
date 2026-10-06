function sigma = bistatic_echo_width_te_pec(phi, k, a)
% truncation order of series
N = 50;
n = -N:1:N;
sigma = zeros(size(phi));

for p = 1:numel(phi)
    sigma(p) = 4*a/k/a * abs(sum((-1j).^(-n) .* dbesselj(n, k*a) .* exp(1j*(2*n+1)*pi/4) ./ dbesselh(n,2,k*a) .* exp(1j*n*phi(p)))).^2;
end
end