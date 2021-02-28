clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')

global C
global CuCond NoCond
global nx ny
global stepX stepY

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²

nx = 75;
ny = 50;

CuCond = 1.7e-8;

%Conductivity map

cMap = zeros(nx,ny);

for i = 1:nx
    for j = 1: ny
        cMap(i,j) = CuCond;
    end
end

G = sparse(nx*ny, nx*ny);
F = zeros(1, nx*ny);

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;
        
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
             F(n) = 1;
            
        elseif i == nx
             G(n,:) = 0;
             G(n,n) = 1;

             F(n) = 1;
          
        
        elseif j == 1
            G(n,:) = 0;
            G(n,n) = 1;
             F(n) = 0;
        elseif j == ny
            G(n,:) = 0;
            G(n,n) = 1;
             F(n) = 0;
        else
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nym = (j-1) + (i-1) * ny;
            nyp = (j+1) + (i-1) * ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm + rxp + rym + ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nym) = rym;
            G(n, nyp) = ryp;
        end
    end
end

% V using Finite Difference Method

V = G\F';

for i = 1:nx
    for j = 1:ny
        n = j + (i-1) * ny;
        VG(i,j) = V(n);
    end   
end

subplot(2,1,1);
set(surf(VG),'linestyle', 'none');

% V using the Analytical Series Method

N =101;

a = ny;
b = nx/2;
VpSig = zeros(nx,ny);

x = linspace(-b,b,nx);
y = linspace(0,a,ny);
% [X,Y] = meshgrid(x,y);

VSig = zeros(nx,ny);
for k = 1:2:N
    for i = 1:nx
        for j = 1:ny
            VpSig(i,j) = ((4*1)/pi) * (1/k) * ((cosh((k*pi*x(i))/a))/(cosh((k*pi*b)/a))) * sin((k*pi*y(j))/a);
        end    
    end
    VSig = (VSig + VpSig);
    subplot(2,1,2);
   
    surf(VSig)
    pause(0.1)
end

