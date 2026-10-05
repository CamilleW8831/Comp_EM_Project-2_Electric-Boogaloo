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

%grid_spacings = [32 20 16 12 8 5 3 2 1 0.5 0.3 0.25]*1e-3; % grid size convergence study
grid_spacings = [1.0]*1e-3;
rmse_sigma_error = zeros(size(grid_spacings));

for gridi = 1:numel(grid_spacings)
dx = grid_spacings(gridi);
dy = dx;
%dx = 2e-3; dy = dx; % spacing [m]
problem_width = 0.08; % problem width excluding PML cells [m]
huygens_distance = 0.02; % distance of Huygens surface from center of grid

PML_width = 0.0150;% width of PML region

PML_cells = ceil(PML_width/dx); % number of cells corrsponding to PML width
PML_atten_dB = 80; % desired PML attenuation (dB)
PML_order = 2; % polynomial order of PML tapering

Nx = floor((problem_width + PML_cells*2*dx) / dx);
Ny = Nx;


Huygens_cells = floor(huygens_distance/dx);
% (i0, j0) is the index of the bottom left corner of the Huygen surface
imin = floor(Nx/2)-Huygens_cells; 
jmin = floor(Ny/2)-Huygens_cells;
% (i1, j1) is the index of the top right corner of the Huygen surface
imax = floor(Nx/2)+Huygens_cells;
jmax = floor(Ny/2)+Huygens_cells;

NFFF_Npts = 256; % number of points for near field to far field transform
NFFF_radius = 3.5e-2; % radius of near field to far field sphere

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

Nt = 1000; % number of time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize the field and material values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Incident fields
Ei_z = zeros(Nx, Ny, Nt); % E-field
Ei_sx = zeros(Nx, Ny, Nt); % split E-field x component (PML)
Ei_sy = zeros(Nx, Ny, Nt); % split E-field y component (PML)

Hi_x = zeros(Nx, Ny, Nt); Hi_y = zeros(Nx, Ny, Nt); % H-field

%%% Total fields
Ez = zeros(Nx, Ny, Nt); % E-field
Ez_sx = zeros(Nx, Ny, Nt); % split E-field x component (PML)
Ez_sy = zeros(Nx, Ny, Nt); % split E-field y component (PML)

Hx = zeros(Nx, Ny, Nt); Hy = zeros(Nx, Ny, Nt); % H-field
% Jz = zeros(Nx, Ny, Nt); % Optional: current

%%% Material values
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
PML_right = X > (x(end)-PML_width);
sigma_x(PML_left) = sigma_max * ((PML_width-X(PML_left)) / PML_width).^PML_order;
sigma_x(PML_right) = sigma_max * ((PML_width-x(end)+X(PML_right)) / PML_width).^PML_order;
PML_top = Y <= PML_width;
PML_bottom = Y > (y(end)-PML_width);
sigma_y(PML_top) = sigma_max * ((PML_width-Y(PML_top)) / PML_width).^PML_order;
sigma_y(PML_bottom) = sigma_max * ((PML_width-y(end)+Y(PML_bottom)) / PML_width).^PML_order;

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

doPEC = true;
% Dielectric circle in the center:

% Coordinates of approximately the center
x0 = x(floor(Nx/2));
y0 = y(floor(Ny/2));

D = 1.5e-2; % diameter [m]
R = D/2; % radius [m]

circle = (X-x0).^2 + (Y-y0).^2 <= R^2;

eps_r(circle) = 9; % assign the permittivity to the circle

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

% source location
ps = imin - 1;  
qs = 1:Ny;

% tapered sinusoid
f0 = 10e9; % frequency
w0 = 2*pi*f0; % angular frequency
sigma_t = 2/f0; % taper time constant

tE = (0:Nt-1)*dt; % time vector, electric field reference

tOffset = Nt/4*dt;

ft = exp(-((tE-tOffset)/sigma_t).^2) .* cos(w0*tE); % excite with modulated Gaussian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build/Update Equations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these values are used in the update for Ez
beta_x = eps0*eps_r/dt + sigma_x/2;
beta_y = eps0*eps_r/dt + sigma_y/2;

alpha_x = eps0*eps_r/dt - sigma_x/2;
alpha_y = eps0*eps_r/dt - sigma_y/2;

% vacuum coefficients for the incident grid
beta0_x = eps0/dt + sigma_x/2;
beta0_y = eps0/dt + sigma_y/2; 

alpha0_x = eps0/dt - sigma_x/2;
alpha0_y = eps0/dt - sigma_y/2;

wb = waitbar(0, 'Running FDTD');
for h = 1:Nt-1 % Nt time steps

    % Update Hx
    for p = 1:Nx
        for q = 1:Ny-1

            % Hx total
            Hx(p,q,h+1) = (alpha_y(p,q) * Hx(p,q,h) - eps_r(p,q)*eps0/(mu0*mu_r(p,q)*dy) * ...
                (Ez(p,q+1,h) - Ez(p,q,h))) / beta_y(p,q);

            % Hx incident (in a vacuum)
            Hi_x(p,q,h+1) = (alpha0_y(p,q) * Hi_x(p,q,h) - eps0/(mu0*dy) * ...
                (Ei_z(p,q+1,h) - Ei_z(p,q,h))) / (beta0_y(p,q));
        end
    end

    % along the bottom edge of the surface whenever q = jmin-1, we will
    % require the term Ez at (p, jmin). This term is corrected such that:
    % Ez_new(p, jmin) = Ez_old(p, jmin) - Ei_z(p, jmin)
    correction_jmin = eps_r(imin:imax,jmin-1).*eps0 ./ (mu0*mu_r(imin:imax,jmin-1)*dy) ...
        .* Ei_z(imin:imax,jmin,h) ./ beta_y(imin:imax,jmin-1);
    %%%
    % apply the jmin correction
    Hx(imin:imax,jmin-1,h+1) = Hx(imin:imax,jmin-1,h+1) + correction_jmin;

    % along the top edge of the surface whenever q = jmax, we will require
    % the term Ez at (p, jmax). This term is corrected such that:
    % Ez_new(p, jmax) = Ez_old(p, jmax) - Ei_z(p, jmax)
    correction_jmax = -eps_r(imin:imax,jmax).*eps0 ./ (mu0*mu_r(imin:imax,jmax)*dy) ...
        .* Ei_z(imin:imax,jmax,h) ./ beta_y(imin:imax,jmax);
    %%%
    % apply the jmax correction
    Hx(imin:imax,jmax,h+1) = Hx(imin:imax,jmax,h+1) + correction_jmax;

    % Update Hy
    for p = 1:Nx-1
        for q = 1:Ny

            % Hy total
            Hy(p,q,h+1) = (alpha_x(p,q) * Hy(p,q,h) + eps_r(p,q)*eps0/(mu0*mu_r(p,q)*dx) * ...
                (Ez(p+1,q,h) - Ez(p,q,h))) / beta_x(p,q);
        
            % Hy incident (in a vacuum)
            Hi_y(p,q,h+1) = (alpha0_x(p,q) * Hi_y(p,q,h) + eps0/(mu0*dx) * ...
                (Ei_z(p+1,q,h) - Ei_z(p,q,h))) / beta0_x(p,q);
        end
    end

    % along the left edge of the surface whenever p = imin-1, we will require
    % the term Ez at (imin, q). This term is corrected such that:
    % Ez_new(imin, q) = Ez_old(imin, q) - Ei_z(imin, q)
    correction_imin = -eps_r(imin-1,jmin:jmax).*eps0 ./ (mu0*mu_r(imin-1,jmin:jmax)*dx) ...
        .* Ei_z(imin,jmin:jmax,h) ./ beta_x(imin-1,jmin:jmax);
    %%%
    % apply the imin correction
    Hy(imin-1,jmin:jmax,h+1) = Hy(imin-1,jmin:jmax,h+1) + correction_imin;

    % along the right edge of the surface whenever p = imax, we will require
    % the term Ez at (imax, q). This term is corrected such that:
    % Ez_new(imax, q) = Ez_old(imax, q) - Ei_z(imax, q)
    correction_imax = eps_r(imax,jmin:jmax).*eps0 ./ (mu0*mu_r(imax,jmin:jmax)*dx) ...
        .* Ei_z(imax,jmin:jmax,h) ./ beta_x(imax,jmin:jmax);
    %%%%
    % apply the imax correction
    Hy(imax,jmin:jmax,h+1) = Hy(imax,jmin:jmax,h+1) + correction_imax;

    % Update Ez
    for p = 2:Nx-1
        for q = 2:Ny-1

            %%%
            % Ez total
            Ez_sx(p,q,h+1) = (alpha_x(p,q) * Ez_sx(p,q,h) ...
                + (Hy(p,q,h+1) - Hy(p-1,q,h+1))/dx) / beta_x(p,q);

            Ez_sy(p,q,h+1) = (alpha_y(p,q) * Ez_sy(p,q,h) ...
                - (Hx(p,q,h+1) - Hx(p,q-1,h+1))/dy) / beta_y(p,q);
        
            % combine Ez total
            Ez(p,q,h+1) = Ez_sx(p,q,h+1) + Ez_sy(p,q,h+1);

            %%%
            % Ez incident (vacuum)
            Ei_sx(p,q,h+1) = (alpha0_x(p,q) * Ei_sx(p,q,h) ...
                + (Hi_y(p,q,h+1) - Hi_y(p-1,q,h+1))/dx) / beta0_x(p,q);

            Ei_sy(p,q,h+1) = (alpha0_y(p,q) * Ei_sy(p,q,h) ...
                - (Hi_x(p,q,h+1) - Hi_x(p,q-1,h+1))/dy) / beta0_y(p,q);

            % combine Ez incident
            Ei_z(p,q,h+1) = Ei_sx(p,q,h+1) + Ei_sy(p,q,h+1);

        end
    end
    
    if doPEC
        Ez(:,:,h+1) = Ez(:,:,h+1) .* ~(circle);
    end

    % along the left edge of the surface whenever p = imin, we will require
    % the term Hy at (imin-1, q) for Ez_sx. This term is corrected such that:
    % Hy_new(imin-1, q) = Hy_old(imin-1, q) + Hi_y(imin-1, q)
    correction_Ez_imin = -Hi_y(imin-1,jmin:jmax,h+1) ./ (dx*beta_x(imin,jmin:jmax));
    %%%
    % apply the imin correction
    Ez_sx(imin,jmin:jmax,h+1) = Ez_sx(imin,jmin:jmax,h+1) + correction_Ez_imin;
    Ez(imin,jmin:jmax,h+1) = Ez(imin,jmin:jmax,h+1) + correction_Ez_imin;

    % along the right edge of the surface whenever p = imax, we will require
    % the term Hy at (imax, q) for Ez_sx. This term is corrected such that:
    % Hy_new(imax, q) = Hy_old(imax, q) + Hi_y(imax, q)
    correction_Ez_imax = Hi_y(imax,jmin:jmax,h+1) ./ (dx*beta_x(imax,jmin:jmax));
    %%%
    % apply the imax correction
    Ez_sx(imax,jmin:jmax,h+1) = Ez_sx(imax,jmin:jmax,h+1) + correction_Ez_imax;
    Ez(imax,jmin:jmax,h+1) = Ez(imax,jmin:jmax,h+1) + correction_Ez_imax;

    % along the bottom edge of the surface whenever q = jmin, we will require
    % the term Hx at (p, jmin-1) for Ez_sy. This term is corrected such that:
    % Hx_new(p, jmin-1) = Hx_old(p, jmin-1) + Hi_x(p, jmin-1)
    correction_Ez_jmin = Hi_x(imin:imax,jmin-1,h+1) ./ (dy*beta_y(imin:imax,jmin));
    %%%
    % apply the jmin correction
    Ez_sy(imin:imax,jmin,h+1) = Ez_sy(imin:imax,jmin,h+1) + correction_Ez_jmin;
    Ez(imin:imax,jmin,h+1) = Ez(imin:imax,jmin,h+1) + correction_Ez_jmin;

    % along the top edge of the surface whenever q = jmax, we will require
    % the term Hx at (p, jmax) for for Ez_sy. This term is corrected such that:
    % Hx_new(p, jmax) = Hx_old(p, jmax) + Hi_x(p, jmax)
    correction_Ez_jmax = -Hi_x(imin:imax,jmax,h+1) ./ (dy*beta_y(imin:imax,jmax));
    %%%
    % apply the jmax correction
    Ez_sy(imin:imax,jmax,h+1) = Ez_sy(imin:imax,jmax,h+1) + correction_Ez_jmax;
    Ez(imin:imax,jmax,h+1) = Ez(imin:imax,jmax,h+1) + correction_Ez_jmax;

    % add the source term (only to the incident field)
    Ei_sx(ps,qs,h+1) =  Ei_sx(ps,qs,h+1) + 0.5*ft(h+1);
    Ei_sx(ps,qs,h+1) =  Ei_sx(ps,qs,h+1) + 0.5*ft(h+1);
    Ei_z(ps,qs,h+1)  = Ei_sx(ps,qs,h+1) + Ei_sy(ps,qs,h+1);

    if mod(h,20) == 0   % update in steps 
        waitbar(h/(Nt-1), wb, sprintf('Iterating: %.1f%%', 100*h/(Nt-1)));
    end

end
close (wb);

%%%%%%%%%%%%%%%%%%
%%% Far Fields %%%
%%%%%%%%%%%%%%%%%%

% compute frequency domain excitation
ft_freq = fftshift(fft(ft));
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
NF_spectrum = fftshift(fft2(Ez_NF)) ./ numel(phi_samp);
% cylindrical harmonic index
N = (-NFFF_Npts/2:NFFF_Npts/2-1);

% propagate cylindrical harmonics to the far field
% this is essentially just the asymptotic expansion of the Hankel function
H = zeros(size(NF_spectrum));
for ni = 1:numel(N)
    n = N(ni);
    H(ni,:) = ((1j)^(-n))*(1./besselh(n, 2, k_F*NFFF_radius)) .* sqrt(k_F) ./ 2;
end
% hankel functions at dc (k=0) give NaN, so we replace with zeros
H(isnan(H)) = 0.0;

FF_spectrum = H .* NF_spectrum;
% convert to spatial far field
FF = ifft(ifftshift(FF_spectrum,1),[],1);
FF = FF ./ abs(ft_freq);

% get index of closest computed frequency to f0
[~, findex] = min(abs(F-f0));
k0 = k_F(findex);
% bistatic echo width (4*pi^2 corrects for FFT normalization)
sigma_sim = 4*pi^2*abs(FF(:,findex)).^2;
sigma_anal = bistatic_echo_width_tm_pec(phi_samp, k0, R);
% compute rmse error in echo width
rmse_sigma = sqrt(sum((sigma_sim.' - sigma_anal).^2)) ./ sqrt(sum(sigma_anal.^2));
rmse_sigma_error(gridi) = rmse_sigma;
end

% figure;
% plot(phi_samp*180/pi, 10*log10(sigma_sim)); hold on;
% plot(phi_samp*180/pi, 10*log10(sigma_ana));
% legend("Simulated", "Analytic");
% xlabel("\phi (deg)");
% ylabel("\sigma(\phi) (dBm)");

%%%%%%%%%%%%%%%%
%%% Frequency Domain Solution %%%
%%%%%%%%%%%%%%%%

tiledlayout(1, 3);

% compute freq domain of field samples
Ez_f = fftshift(fft(Ez, [], 3), 3);
Ezi_f = fftshift(fft(Ei_z, [], 3), 3);

% incident field
nexttile;
imagesc(x, y, real(Ezi_f(:,:,findex)')); hold on;
colormap magma;
cl = 0.50*max(abs(real(Ezi_f(:,:,findex)')), [], 'all'); % more extreme color gradient for clear visual
clim([-cl cl]);
axis equal tight;
title("Incident E_z");
plot(x0 + R*cos(linspace(0,2*pi,300)), y0 + R*sin(linspace(0,2*pi,300)), 'k', 'LineWidth', 1.5);

plot(x_samp, y_samp, "g:");
% Optional: Huygens surface
rectangle('Position', [x(imin), y(jmin), x(imax)-x(imin), y(jmax)-y(jmin)], ...
    'EdgeColor', 'w', 'LineStyle', '--', 'LineWidth', 1.5);
% Optional: PML inner boundary
rectangle('Position', [x(PML_cells+1), y(PML_cells+1), ...
    x(Nx-PML_cells)-x(PML_cells+1), y(Ny-PML_cells)-y(PML_cells+1)], ...
    'EdgeColor', 'w', 'LineStyle', ':', 'LineWidth', 1.5);

% scattered field
mask = (X/dx >= imin-1 & X/dx < imax & Y/dy >= jmin-1 & Y/dy < jmax);
Ez_f = Ez_f - (mask.*Ezi_f);

nexttile;
imagesc(x, y, real(Ez_f(:,:,findex)')); hold on;
colormap magma;
cl = 0.50*max(abs(real(Ez_f(:,:,findex)')), [], 'all'); % more extreme color gradient for clear visual
clim([-cl cl]);
axis equal tight;
title("Scattered E_z");
plot(x0 + R*cos(linspace(0,2*pi,300)), y0 + R*sin(linspace(0,2*pi,300)), 'k', 'LineWidth', 1.5);

plot(x_samp, y_samp, "g:");
% Optional: Huygens surface
rectangle('Position', [x(imin), y(jmin), x(imax)-x(imin), y(jmax)-y(jmin)], ...
    'EdgeColor', 'w', 'LineStyle', '--', 'LineWidth', 1.5);
% Optional: PML inner boundary
rectangle('Position', [x(PML_cells+1), y(PML_cells+1), ...
    x(Nx-PML_cells)-x(PML_cells+1), y(Ny-PML_cells)-y(PML_cells+1)], ...
    'EdgeColor', 'w', 'LineStyle', ':', 'LineWidth', 1.5);

% total field
Ez_t = Ez_f + Ezi_f;
nexttile;
imagesc(x, y, real(Ez_t(:,:,findex)')); hold on;
c = colorbar;
ylabel(c, "Re(E_z)");
colormap magma;
cl = 0.50*max(abs(real(Ez_t(:,:,findex)')), [], 'all'); % more extreme color gradient for clear visual
clim([-cl cl]);
axis equal tight;
title("Total E_z");
plot(x0 + R*cos(linspace(0,2*pi,300)), y0 + R*sin(linspace(0,2*pi,300)), 'k', 'LineWidth', 1.5);

plot(x_samp, y_samp, "g:");
% Optional: Huygens surface
rectangle('Position', [x(imin), y(jmin), x(imax)-x(imin), y(jmax)-y(jmin)], ...
    'EdgeColor', 'w', 'LineStyle', '--', 'LineWidth', 1.5);
% Optional: PML inner boundary
rectangle('Position', [x(PML_cells+1), y(PML_cells+1), ...
    x(Nx-PML_cells)-x(PML_cells+1), y(Ny-PML_cells)-y(PML_cells+1)], ...
    'EdgeColor', 'w', 'LineStyle', ':', 'LineWidth', 1.5);


%%%%%%%%%%%%%%%%
%%% Playback %%%
%%%%%%%%%%%%%%%%

figure();

for h = 1:6:Nt

    % plot
    imagesc(x, y, Ez(:,:,h)'); hold on;
    plot(x0 + R*cos(linspace(0,2*pi,300)), y0 + R*sin(linspace(0,2*pi,300)), 'k', 'LineWidth', 1.5);

    plot(x_samp, y_samp, "g:");

    % Optional: Huygens surface
    rectangle('Position', [x(imin), y(jmin), x(imax)-x(imin), y(jmax)-y(jmin)], ...
        'EdgeColor', 'w', 'LineStyle', '--', 'LineWidth', 1.5);

    % Optional: PML inner boundary
    rectangle('Position', [x(PML_cells+1), y(PML_cells+1), ...
        x(Nx-PML_cells)-x(PML_cells+1), y(Ny-PML_cells)-y(PML_cells+1)], ...
        'EdgeColor', 'w', 'LineStyle', ':', 'LineWidth', 1.5);

    % % Optional: source location
    % plot(x(ps), y(qs), 'wp', 'MarkerFaceColor', 'k', 'MarkerSize', 12);

    % colorbar settings
    colorbar;
    colormap magma;
    cl = 0.50*max(abs(Ez(:))); % more extreme color gradient for clear visual
    clim([-cl cl]);

    % axis & figure settings
    set(gca,'YDir','normal', 'TickLength', [0,0], 'FontName', 'Times', 'FontSize', 18);
    set(gcf, 'Color', 'w')
    axis image;

    % labeling
    title("E [V/m] at t = " + round((h-1)*dt*1e9) + " ns"); 
    xlabel("x [m]")
    ylabel("y [m]")

    drawnow; % animate

end
