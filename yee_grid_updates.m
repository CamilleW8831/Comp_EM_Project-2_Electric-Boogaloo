clc; clear; close all;

%%%%%%%%%%%%%%%%%
%%% Constants %%%
%%%%%%%%%%%%%%%%%

c = 3* 1e8; % speed of light in vacuum [m/s]
eps0 = 8.85*1e-12; % permativitty of free space [F/m]
mu0 = 4*pi*1e-7; % permeability of free space [H/m]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spatial Indexing %%%
%%%%%%%%%%%%%%%%%%%%%%%%

Nx = 101; Ny = 101; % number of points
dx = 1; dy = 1; % spacing

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

% time step s.t. stability criterion is met
dt = 0.9 / (c * sqrt(1/dx^2 + 1/dy^2));

Nt = 500; % number of time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize the field and material values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ez = zeros(Nx, Ny, Nt); % E-field
Hx = zeros(Nx, Ny, Nt); Hy = zeros(Nx, Ny, Nt); % H-field
Jz = zeros(Nx, Ny, Nt); % Current

eps_r = ones(Nx, Ny); % permittivity 
mu_r = ones(Nx, Ny); % permeability
sigma = zeros(Nx, Ny); % conductivity

%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Geometry %%%
%%%%%%%%%%%%%%%%%%%%%%%

% Dielectric cirlce in the center:

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

ps  = round(Nx/4);  qs = round(Ny/2); % source loacation

% tapered sinusoid
f0 = 10e6; % frequency
w0 = 2*pi*f0; % angular frequency
sigma_t = 1/f0; % taper time constant

tE = (0:Nt-1)*dt; % time vector, electric field reference
ft = (1 - exp(-tE/sigma_t)) .* sin(w0*tE); % exciting function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build/Update Equations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these values are used in the update for Ez
beta  = eps0*eps_r/dt + sigma/2;
alpha = eps0*eps_r/dt - sigma/2;

% Extremity field values are not updated from 0, implying PEC boundary conditions

for h = 1:Nt-1 % Nt time steps

    % Update Hx
    for p = 1:Nx
        for q = 1:Ny-1
            Hx(p,q,h+1) = Hx(p,q,h) - dt/(mu0*mu_r(p,q)*dy) * (Ez(p,q+1,h) - Ez(p,q,h));
        end
    end

    % Update Hy
    for p = 1:Nx-1
        for q = 1:Ny
            Hy(p,q,h+1) = Hy(p,q,h) + dt/(mu0*mu_r(p,q)*dx) * (Ez(p+1,q,h) - Ez(p,q,h));
        end
    end

    % Update Ez
    for p = 2:Nx-1
        for q = 2:Ny-1
            curlH = (Hy(p,q,h+1) - Hy(p-1,q,h+1))/dx ...
                  - (Hx(p,q,h+1) - Hx(p,q-1,h+1))/dy;
            Ez(p,q,h+1) = (alpha(p,q) * Ez(p,q,h) + curlH - Jz(p,q,h+1)) / beta(p,q);
        end
    end

    % add the source term
    Ez(ps,qs,h+1) = Ez(ps,qs,h+1) + ft(h+1);

end

%%%%%%%%%%%%%%%%
%%% Playback %%%
%%%%%%%%%%%%%%%%

figure();

for h = 1:5:Nt % every 5th time step

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






