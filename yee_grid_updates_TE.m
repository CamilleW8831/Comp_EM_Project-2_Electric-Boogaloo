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

PML_cells = 45; % number of cells corrsponding to PML width
PML_width = PML_cells*dx; % width of PML region
PML_atten_dB = 80; % desired PML attenuation (dB)
PML_order = 2; % polynomial order of PML tapering

PML_Huygen_cells = 10; % separation between the Huygen surface and the PML [cells]
% (i0, j0) is the index of the bottom left corner of the Huygen surface
imin = PML_cells + PML_Huygen_cells; 
jmin = PML_cells + PML_Huygen_cells;
% (i1, j1) is the index of the top right corner of the Huygen surface
imax = Nx - PML_cells - PML_Huygen_cells;
jmax = Ny - PML_cells - PML_Huygen_cells;

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

%%% Incident fields
Hi_z = zeros(Nx, Ny, Nt); % H-field
Hi_sx = zeros(Nx, Ny, Nt); % split H-field x component (PML)
Hi_sy = zeros(Nx, Ny, Nt); % split H-field y component (PML)

Ei_x = zeros(Nx, Ny, Nt); Ei_y = zeros(Nx, Ny, Nt); % E-field

%%% Total fields
Hz = zeros(Nx, Ny, Nt); % H-field
Hz_sx = zeros(Nx, Ny, Nt); % split H-field x component (PML)
Hz_sy = zeros(Nx, Ny, Nt); % split H-field y component (PML)

Ex = zeros(Nx, Ny, Nt); Ey = zeros(Nx, Ny, Nt); % H-field
% Jz = zeros(Nx, Ny, Nt); % Optional: current

%%% Material values
eps_r = ones(Nx, Ny); % permittivity 
mu_r = ones(Nx, Ny); % permeability

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create PML Layers %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% compute maximum loss for desired attenuation
sigma_max = -(PML_order+1)/(2*eta0*PML_width)*log(10^(-PML_atten_dB/20)) * mu0/eps0;

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

% Dielectric circle in the center:

% Coordinates of approximately the center
x0 = floor(x(end)/2);
y0 = floor(y(end)/2); 

D = 20; % diameter [m]
R = D/2; % radius [m]

dielectric = (X-x0).^2 + (Y-y0).^2 <= R^2;
% eps_r(dielectric) = 9; % assign the permittivity to the circle
Ex_PEC = (X-x0).^2 + (Y+dy/2-y0).^2 <= R^2;
Ey_PEC = (X+dx/2-x0).^2 + (Y-y0).^2 <= R^2;

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
ps = imin - round(PML_Huygen_cells/2);  
qs = 1:Ny;

% tapered sinusoid
f0 = 20e6; % frequency
w0 = 2*pi*f0; % angular frequency
sigma_t = 1/f0; % taper time constant

tE = (0:Nt-1)*dt; % time vector, electric field reference
ft = (1 - exp(-tE/sigma_t)) .* sin(w0*tE); % exciting function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build/Update Equations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these values are used in the update for Hz
beta_x = mu0*mu_r/dt + sigma_x/2;
beta_y = mu0*mu_r/dt + sigma_y/2;

alpha_x = mu0*mu_r/dt - sigma_x/2;
alpha_y = mu0*mu_r/dt - sigma_y/2;

% vacuum coefficients for the incident grid
beta0_x = mu0/dt + sigma_x/2;
beta0_y = mu0/dt + sigma_y/2; 

alpha0_x = mu0/dt - sigma_x/2;
alpha0_y = mu0/dt - sigma_y/2;

wb = waitbar(0, 'Running FDTD');
for h = 1:Nt-1 % Nt time steps

    % Update Ex
    for p = 1:Nx
        for q = 1:Ny-1

            % Ex total
            Ex(p,q,h+1) = (alpha_y(p,q) * Ex(p,q,h) + mu_r(p,q)*mu0/(eps0*eps_r(p,q)*dy) * ...
                (Hz(p,q+1,h) - Hz(p,q,h))) / beta_y(p,q);

            % Ex incident (in a vacuum)
            Ei_x(p,q,h+1) = (alpha0_y(p,q) * Ei_x(p,q,h) + mu0/(eps0*dy) * ...
                (Hi_z(p,q+1,h) - Hi_z(p,q,h))) / (beta0_y(p,q));
        end
    end

    % along the bottom edge of the surface whenever q = jmin-1, we will
    % require the term Hz at (p, jmin). This term is corrected such that:
    % Hz_new(p, jmin) = Hz_old(p, jmin) - Hi_z(p, jmin)
    correction_jmin = -mu_r(imin:imax,jmin-1).*mu0 ./ (eps0*eps_r(imin:imax,jmin-1)*dy) ...
        .* Hi_z(imin:imax,jmin,h) ./ beta_y(imin:imax,jmin-1);
    %%%
    % apply the jmin correction
    Ex(imin:imax,jmin-1,h+1) = Ex(imin:imax,jmin-1,h+1) + correction_jmin;

    % along the top edge of the surface whenever q = jmax, we will require
    % the term Hz at (p, jmax). This term is corrected such that:
    % Hz_new(p, jmax) = Hz_old(p, jmax) - Hi_z(p, jmax)
    correction_jmax = mu_r(imin:imax,jmax).*mu0 ./ (eps0*eps_r(imin:imax,jmax)*dy) ...
        .* Hi_z(imin:imax,jmax,h) ./ beta_y(imin:imax,jmax);
    %%%
    % apply the jmax correction
    Ex(imin:imax,jmax,h+1) = Ex(imin:imax,jmax,h+1) + correction_jmax;

    % Update Ey
    for p = 1:Nx-1
        for q = 1:Ny

            % Ey total
            Ey(p,q,h+1) = (alpha_x(p,q) * Ey(p,q,h) - mu_r(p,q)*mu0/(eps0*eps_r(p,q)*dx) * ...
                (Hz(p+1,q,h) - Hz(p,q,h))) / beta_x(p,q);
        
            % Ey incident (in a vacuum)
            Ei_y(p,q,h+1) = (alpha0_x(p,q) * Ei_y(p,q,h) - mu0/(eps0*dx) * ...
                (Hi_z(p+1,q,h) - Hi_z(p,q,h))) / beta0_x(p,q);
        end
    end

    % along the left edge of the surface whenever p = imin-1, we will require
    % the term Hz at (imin, q). This term is corrected such that:
    % Hz_new(imin, q) = Hz_old(imin, q) - Hi_z(imin, q)
    correction_imin = mu_r(imin-1,jmin:jmax).*mu0 ./ (eps0*eps_r(imin-1,jmin:jmax)*dx) ...
        .* Hi_z(imin,jmin:jmax,h) ./ beta_x(imin-1,jmin:jmax);
    %%%
    % apply the imin correction
    Ey(imin-1,jmin:jmax,h+1) = Ey(imin-1,jmin:jmax,h+1) + correction_imin;

    % along the right edge of the surface whenever p = imax, we will require
    % the term Hz at (imax, q). This term is corrected such that:
    % Hz_new(imax, q) = Hz_old(imax, q) - Hi_z(imax, q)
    correction_imax = -mu_r(imax,jmin:jmax).*mu0 ./ (eps0*eps_r(imax,jmin:jmax)*dx) ...
        .* Hi_z(imax,jmin:jmax,h) ./ beta_x(imax,jmin:jmax);
    %%%%
    % apply the imax correction
    Ey(imax,jmin:jmax,h+1) = Ey(imax,jmin:jmax,h+1) + correction_imax;

    % select the current time step of Ex and Ey
    Ex_tan = Ex(:,:,h+1); Ey_tan = Ey(:,:,h+1);

    % Set the values inside of the PEC = 0
    Ex_tan(Ex_PEC) = 0; Ex(:,:,h+1) = Ex_tan;
    Ey_tan(Ey_PEC) = 0;  Ey(:,:,h+1) = Ey_tan;

    % Update Hz
    for p = 2:Nx-1
        for q = 2:Ny-1

            %%%
            % Hz total
            Hz_sx(p,q,h+1) = (alpha_x(p,q) * Hz_sx(p,q,h) ...
                + (Ey(p-1,q,h+1) - Ey(p,q,h+1))/dx) / beta_x(p,q);

            Hz_sy(p,q,h+1) = (alpha_y(p,q) * Hz_sy(p,q,h) ...
                - (Ex(p,q-1,h+1) - Ex(p,q,h+1))/dy) / beta_y(p,q);
        
            % combine Hz total
            Hz(p,q,h+1) = Hz_sx(p,q,h+1) + Hz_sy(p,q,h+1);

            %%%
            % Hz incident (vacuum)
            Hi_sx(p,q,h+1) = (alpha0_x(p,q) * Hi_sx(p,q,h) ...
                + (Ei_y(p-1,q,h+1) - Ei_y(p,q,h+1))/dx) / beta0_x(p,q);

            Hi_sy(p,q,h+1) = (alpha0_y(p,q) * Hi_sy(p,q,h) ...
                - (Ei_x(p,q-1,h+1) - Ei_x(p,q,h+1))/dy) / beta0_y(p,q);

            % combine Hz incident
            Hi_z(p,q,h+1) = Hi_sx(p,q,h+1) + Hi_sy(p,q,h+1);

        end
    end

    % along the left edge of the surface whenever p = imin, we will require
    % the term Ey at (imin-1, q) for Hz_sx. This term is corrected such that:
    % Ey_new(imin-1, q) = Ey_old(imin-1, q) + Ei_y(imin-1, q)
    correction_Hz_imin = Ei_y(imin-1,jmin:jmax,h+1) ./ (dx*beta_x(imin,jmin:jmax));
    %%%
    % apply the imin correction
    Hz_sx(imin,jmin:jmax,h+1) = Hz_sx(imin,jmin:jmax,h+1) + correction_Hz_imin;
    Hz(imin,jmin:jmax,h+1) = Hz(imin,jmin:jmax,h+1) + correction_Hz_imin;

    % along the right edge of the surface whenever p = imax, we will require
    % the term Ey at (imax, q) for Hz_sx. This term is corrected such that:
    % Ey_new(imax, q) = Ey_old(imax, q) + Ei_y(imax, q)
    correction_Hz_imax = -Ei_y(imax,jmin:jmax,h+1) ./ (dx*beta_x(imax,jmin:jmax));
    %%%
    % apply the imax correction
    Hz_sx(imax,jmin:jmax,h+1) = Hz_sx(imax,jmin:jmax,h+1) + correction_Hz_imax;
    Hz(imax,jmin:jmax,h+1) = Hz(imax,jmin:jmax,h+1) + correction_Hz_imax;

    % along the bottom edge of the surface whenever q = jmin, we will require
    % the term Ex at (p, jmin-1) for Hz_sy. This term is corrected such that:
    % Ex_new(p, jmin-1) = Ex_old(p, jmin-1) + Ei_x(p, jmin-1)
    correction_Hz_jmin = -Ei_x(imin:imax,jmin-1,h+1) ./ (dy*beta_y(imin:imax,jmin));
    %%%
    % apply the jmin correction
    Hz_sy(imin:imax,jmin,h+1) = Hz_sy(imin:imax,jmin,h+1) + correction_Hz_jmin;
    Hz(imin:imax,jmin,h+1) = Hz(imin:imax,jmin,h+1) + correction_Hz_jmin;

    % along the top edge of the surface whenever q = jmax, we will require
    % the term Ex at (p, jmax) for for Hz_sy. This term is corrected such that:
    % Ex_new(p, jmax) = Ex_old(p, jmax) + Ei_x(p, jmax)
    correction_Hz_jmax = Ei_x(imin:imax,jmax,h+1) ./ (dy*beta_y(imin:imax,jmax));
    %%%
    % apply the jmax correction
    Hz_sy(imin:imax,jmax,h+1) = Hz_sy(imin:imax,jmax,h+1) + correction_Hz_jmax;
    Hz(imin:imax,jmax,h+1) = Hz(imin:imax,jmax,h+1) + correction_Hz_jmax;

    % add the source term (only to the incident field, split evenly)
    Hi_sx(ps,qs,h+1) = Hi_sx(ps,qs,h+1) + ft(h+1);
    Hi_z(ps,qs,h+1)  = Hi_sx(ps,qs,h+1) + Hi_sy(ps,qs,h+1);

    if mod(h,20) == 0   % update in steps 
        waitbar(h/(Nt-1), wb, sprintf('Iterating: %.1f%%', 100*h/(Nt-1)));
    end

end
close (wb);

% % % % %%%%%%%%%%%%%%%%%%
% % % % %%% Far Fields %%%
% % % % %%%%%%%%%%%%%%%%%%
% % % % 
% % % % % compute frequency domain excitation
% % % % ft_freq = fft(ft);
% % % % % frequency index
% % % % F = (1/dt)/Nt*(-Nt/2:Nt/2-1);
% % % % % wavenumber
% % % % k_F = 2*pi*F/c;
% % % % 
% % % % % compute circular sampling surface
% % % % phi_samp = linspace(0, 2*pi, NFFF_Npts);
% % % % x_samp = NFFF_radius*cos(phi_samp)+x0;
% % % % y_samp = NFFF_radius*sin(phi_samp)+y0;
% % % % % construct appropriate arrays for interpolation over time simultaneously
% % % % time_i = 1:size(Ez,3);
% % % % time_grid = reshape(time_i, 1, 1, Nt);
% % % % time_shape = ones(1, 1, Nt);
% % % % Xn = X .* time_shape;
% % % % Yn = Y .* time_shape;
% % % % Tn = ones(size(Ez,1), size(Ez,2)) .* time_grid;
% % % % Xs = x_samp .* ones(size(Ez,3),1);
% % % % Ys = y_samp .* ones(size(Ez,3),1);
% % % % Ts = time_i.' .* ones(1, NFFF_Npts);
% % % % 
% % % % Ez_NF = interpn(Xn, Yn, Tn, Ez, Xs, Ys, Ts).';
% % % % % convert sampled near fields to cylindrical harmonics and frequency domain
% % % % % both conversions are done simultaneously through 2d fft
% % % % NF_spectrum = fftshift(fft2(Ez_NF));
% % % % % cylindrical harmonic index
% % % % N = (-NFFF_Npts/2:NFFF_Npts/2-1);
% % % % 
% % % % % propagate cylindrical harmonics to the far field
% % % % % this is essentially just the asymptotic expansion of the Hankel function
% % % % H = zeros(size(NF_spectrum));
% % % % for ni = 1:numel(N)
% % % %     n = N(ni);
% % % %     H(ni,:) = (1j^n)*(1./besselh(n, 2, k_F*NFFF_radius));
% % % % end
% % % % % hankel functions at dc (k=0) give NaN, so we replace with zeros
% % % % H(isnan(H)) = 0.0;
% % % % 
% % % % FF_spectrum = H .* NF_spectrum;
% % % % % convert to spatial far field
% % % % FF = ifft2(ifftshift(FF_spectrum));
% % % % FF = FF ./ ft_freq;

%%%%%%%%%%%%%%%%
%%% Playback %%%
%%%%%%%%%%%%%%%%

figure();

for h = 1:3:Nt

    % plot
    imagesc(x, y, Hz(:,:,h)'); hold on;
    plot(x0 + R*cos(linspace(0,2*pi,300)), y0 + R*sin(linspace(0,2*pi,300)), 'k', 'LineWidth', 1.5);

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
    cl = 0.5*max(abs(Hz(:))); % more extreme color gradient for clear visual
    clim([-cl cl]);

    % axis & figure settings
    set(gca,'YDir','normal', 'TickLength', [0,0], 'FontName', 'Times', 'FontSize', 18);
    set(gcf, 'Color', 'w')
    axis image;
    xticks(0:50:Nx);
    yticks(0:50:Ny);

    % labeling
    title("H [A/m] at t = " + round((h-1)*dt*1e9) + " ns"); 
    xlabel("x [m]")
    ylabel("y [m]")

    drawnow; % animate

end