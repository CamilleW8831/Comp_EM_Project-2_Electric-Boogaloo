clc; clear; close all;

%%%%%%%%%%%%%%%%%
%%% Constants %%%
%%%%%%%%%%%%%%%%%

c = 3* 1e8; % speed of light in vacuum [m/s]
eps0 = 8.85*1e-12; % permittivity of free space [F/m]
mu0 = 4*pi*1e-7; % permeability of free space [H/m]
eta0 = sqrt(mu0/eps0); % impedance of free space [ohms]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spatial Indexing %%%
%%%%%%%%%%%%%%%%%%%%%%%%

Nx = 201; Ny = 201; % number of points
dx = 0.5; dy = 0.5; % spacing

PML_width = 45*dx; % width of PML region
PML_atten_dB = 80; % desired PML attenuation (dB)
PML_order = 2; % polynomial order of PML tapering

NFFF_Npts = 256; % number of points for near field to far field transform
NFFF_radius = 30.0; % radius of near field to far field sphere

x = (0:(Nx-1)) * dx; % horizontal index, corresponding to the right corner
y = (0:(Ny-1)) * dy; % vertical index, corresponding to the bottom corner

Nx = size(x,2); % number of points in x
Ny = size(y,2); % number of points in y

[X,Y] = ndgrid(x,y); % form the 2D grid

% % Optional: visualize the grid
% figure(); hold on;
% plot(X,  Y, 'k');
% plot(X', Y', 'k');
% xlabel('x'); 
% ylabel('y');
% title('2D Spatial Grid')
% axis equal tight;
% set(gcf,'Color','w');

%%%%%%%%%%%%%%%%%%%%%
%%% Time indexing %%%
%%%%%%%%%%%%%%%%%%%%%

stability_factor = 0.9; % scale factor on CFL condition (<=1 is stable)
% time step s.t. stability criterion is met
dt = stability_factor / (c * sqrt(1/dx^2 + 1/dy^2));

Nt = 500; % number of time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize the field and material values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ez = zeros(Nx, Ny, Nt); % E-field
Ez_sx = zeros(Nx, Ny, Nt); % split E-field x component (PML)
Ez_sy = zeros(Nx, Ny, Nt); % split E-field y component (PML)

Hx = zeros(Nx, Ny, Nt); Hy = zeros(Nx, Ny, Nt); % H-field
Jz = zeros(Nx, Ny, Nt); % Current

eps_r = ones(Nx, Ny); % permittivity 
mu_r = ones(Nx, Ny); % permeability

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create PML Layers %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% compute maximum loss for desired attenuation
sigma_max = -(PML_order+1)/(2*eta0*PML_width)*log(10^(-PML_atten_dB/20));

% loss for x and y directed fields
sigma_x = zeros(Nx, Ny);
sigma_y = zeros(Nx, Ny);
% construct the four boundary regions which make up the PML
PML_left = X <= PML_width;
PML_right = X > (Nx*dx-PML_width);
sigma_x(PML_left) = sigma_max * ((PML_width-X(PML_left)) / PML_width).^PML_order;
sigma_x(PML_right) = sigma_max * ((PML_width-Nx*dx+X(PML_right)) / PML_width).^PML_order;
PML_top = Y <= PML_width;
PML_bottom = Y > (Ny*dy-PML_width);
sigma_y(PML_top) = sigma_max * ((PML_width-Y(PML_top)) / PML_width).^PML_order;
sigma_y(PML_bottom) = sigma_max * ((PML_width-Ny*dy+Y(PML_bottom)) / PML_width).^PML_order;

% % Optional: check the PML parameters
% figure();
% imagesc(x, y, sigma_x');
% xlabel('x [m]');
% ylabel('y [m]');
% title('\sigma_x');
% axis equal tight;
% set(gca,'YDir','normal');
% set(gcf,'Color','w');
% colorbar;
% 
% figure();
% imagesc(x, y, sigma_y');
% xlabel('x [m]');
% ylabel('y [m]');
% title('\sigma_y');
% axis equal tight;
% set(gca,'YDir','normal');
% set(gcf,'Color','w');
% colorbar;

%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Geometry %%%
%%%%%%%%%%%%%%%%%%%%%%%

% Dielectric circle in the center:

% Coordinates of approximately the center
x0 = floor(x(end)/2);
y0 = floor(y(end)/2); 

D = 20; % diameter [m]
R = D/2; % radius [m]

dielectric = (X-x0).^2 + (Y-y0).^2 <= R^2;
eps_r(dielectric) = 9; % assign the permittivity to the circle

% % Optional: check the geometry
% figure();
% imagesc(x, y, eps_r');
% xlabel('x [m]');
% ylabel('y [m]');
% title('\epsilon_r');
% axis equal tight;
% set(gca,'YDir','normal');
% set(gcf,'Color','w');
% colorbar;
% 
% figure();
% imagesc(x, y, mu_r');
% xlabel('x [m]');
% ylabel('y [m]');
% title('\mu_r');
% axis equal tight;
% set(gca,'YDir','normal');
% set(gcf,'Color','w');
% colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Excite the Problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

ps  = round(Nx/4);  qs = round(Ny/2); % source location

% tapered sinusoid
f0 = 20e6; % frequency
w0 = 2*pi*f0; % angular frequency
sigma_t = 1/f0; % taper time constant

tE = (0:Nt-1)*dt; % time vector, electric field reference
ft = (1 - exp(-tE/sigma_t)) .* sin(w0*tE); % exciting function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build/Update Equations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these values are used in the update for Ez
beta_x = eps0*eps_r/dt + sigma_x/2;
beta_y = eps0*eps_r/dt + sigma_y/2;

alpha_x = eps0*eps_r/dt - sigma_x/2;
alpha_y = eps0*eps_r/dt - sigma_y/2;

% Extremity field values are not updated from 0, implying PEC boundary conditions

for h = 1:Nt-1 % Nt time steps
    % Update Hx
    for p = 1:Nx
        for q = 1:Ny-1
            Hx(p,q,h+1) = (alpha_x(p,q)*Hx(p,q,h) - eps_r(p,q)*eps0/(mu0*mu_r(p,q)*dy) * (Ez(p,q+1,h) - Ez(p,q,h))) / beta_x(p,q);
        end
    end

    % Update Hy
    for p = 1:Nx-1
        for q = 1:Ny
            Hy(p,q,h+1) = (alpha_y(p,q)*Hy(p,q,h) + eps_r(p,q)*eps0/(mu0*mu_r(p,q)*dx) * (Ez(p+1,q,h) - Ez(p,q,h))) / beta_y(p,q);
        end
    end

    % Update Ez
    for p = 2:Nx-1
        for q = 2:Ny-1
            Ez_sx(p,q,h+1) = (alpha_x(p,q) * Ez_sx(p,q,h) + (Hy(p,q,h+1) - Hy(p-1,q,h+1))/dx) / beta_x(p,q);
            Ez_sy(p,q,h+1) = (alpha_y(p,q) * Ez_sy(p,q,h) - (Hx(p,q,h+1) - Hx(p,q-1,h+1))/dy) / beta_y(p,q);
        
            % combine split fields
            Ez(p,q,h+1) = Ez_sx(p,q,h+1) + Ez_sy(p,q,h+1);
        end
    end

    % add the source term
    Ez(ps,qs,h+1) = Ez(ps,qs,h+1) + ft(h+1);

end

%%%%%%%%%%%%%%%%%%
%%% Far Fields %%%
%%%%%%%%%%%%%%%%%%

% compute frequency domain excitation
ft_freq = fft(ft);
% frequency index
F = (1/dt)/Nt*(-Nt/2:Nt/2-1);
% wavenumber
k_F = 2*pi*F/c;

% compute circular sampling surface
phi_samp = linspace(0, 2*pi, NFFF_Npts);
x_samp = NFFF_radius*cos(phi_samp)+x0;
y_samp = NFFF_radius*sin(phi_samp)+y0;
% construct appropriate arrays for interpolation over time simultaneously
time_i = 1:size(Ez,3);
time_grid = reshape(time_i, 1, 1, Nt);
time_shape = ones(1, 1, Nt);
Xn = X .* time_shape;
Yn = Y .* time_shape;
Tn = ones(size(Ez,1), size(Ez,2)) .* time_grid;
Xs = x_samp .* ones(size(Ez,3),1);
Ys = y_samp .* ones(size(Ez,3),1);
Ts = time_i.' .* ones(1, NFFF_Npts);

Ez_NF = interpn(Xn, Yn, Tn, Ez, Xs, Ys, Ts).';
% convert sampled near fields to cylindrical harmonics and frequency domain
% both conversions are done simultaneously through 2d fft
NF_spectrum = fftshift(fft2(Ez_NF));
% cylindrical harmonic index
N = (-NFFF_Npts/2:NFFF_Npts/2-1);

% propagate cylindrical harmonics to the far field
% this is essentially just the asymptotic expansion of the Hankel function
H = zeros(size(NF_spectrum));
for ni = 1:numel(N)
    n = N(ni);
    H(ni,:) = (1j^n)*(1./besselh(n, 2, k_F*NFFF_radius));
end
% hankel functions at dc (k=0) give NaN, so we replace with zeros
H(isnan(H)) = 0.0;

FF_spectrum = H .* NF_spectrum;
% convert to spatial far field
FF = ifft2(ifftshift(FF_spectrum));
FF = FF ./ ft_freq;


%%%%%%%%%%%%%%%%
%%% Playback %%%
%%%%%%%%%%%%%%%%

figure();

for h = 1:3:Nt % every 5th time step

    % plot
    imagesc(x, y, Ez(:,:,h)'); hold on;
    plot(x0 + R*cos(linspace(0,2*pi,300)), y0 + R*sin(linspace(0,2*pi,300)), 'k', 'LineWidth', 1.5);

    % colorbar settings
    colorbar;
    colormap magma;
    cl = 0.25*max(abs(Ez(:))); % more extreme color gradient for clear visual
    clim([-cl cl]);

    % axis & figure settings
    set(gca,'YDir','normal', 'TickLength', [0,0], 'FontName', 'Times', 'FontSize', 18);
    set(gcf, 'Color', 'w')
    axis image;
    xticks(0:20:Nx);
    yticks(0:20:Ny);

    % labelling
    title("t = " + (h-1)*dt*1e9 + " ns"); 
    xlabel("x [m]")
    ylabel("y [m]")

    drawnow; % animate

end






