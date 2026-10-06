function sigma = bistatic_echo_width_te_eps(phi, k, a, ep_r)
% truncation order of series
N = 50;
n = -N:1:N;
sigma = zeros(size(phi));

kd = sqrt(ep_r)*k;

for p = 1:numel(phi)
    alpha_n = (-1j).^(-n) .* (sqrt(ep_r).*dbesselj(n, k*a).*besselj(n, kd*a) - besselj(n, k*a).*dbesselj(n, kd*a)) ./ (sqrt(ep_r).*dbesselh(n,2,k*a).*besselj(n,kd*a) - besselh(n,2,k*a).*dbesselj(n,kd*a));
    sigma(p) = 4*a/k/a * abs(sum(alpha_n .* exp(1j*(2*n+1)*pi/4).* exp(1j*n*phi(p)))).^2;
end
end