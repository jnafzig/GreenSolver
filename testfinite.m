clear;

Nelem = 201;
Nstencil = 7;

x = linspace(-10,10,Nelem)';
dx = x(2)-x(1);

R = 3;
L = 4;
v0 = -1;

% v = 0*x;
% v(x>=(R-L/2)&x<=(R+L/2)) = v0;
v1 = v0*cosh(x-R/2).^(-2);
v2 = v0*cosh(x+R/2).^(-2);

[dn,W,G] = green(Nelem,Nstencil,dx);

E1 = -1;
E2 = -.75;
eta = 1;

% Npath = 30;
% 
% [ZR,ZI] = ndgrid(linspace(E1,E2,Npath),linspace(eta2,eta1,Npath));
% Z = ZR+ZI*1i;

dn1 = @(E) dn(E,v1);
G1 = @(E) G(E,v1);

dn2 = @(E) dn(E,v2);
G2 = @(E) G(E,v2);

dn = @(E) dn(E,v1+v2);
GE = @(E) G(E,v1+v2);

n1 = real(2*1i/pi*integral(dn1,E1,E2,...
    'Waypoints',[E1+eta*1i,E2+eta*1i],...
    'ArrayValued',true));
n2 = real(2*1i/pi*integral(dn2,E1,E2,...
    'Waypoints',[E1+eta*1i,E2+eta*1i],...
    'ArrayValued',true));
n = real(2*1i/pi*integral(dn,E1,E2,...
    'Waypoints',[E1+eta*1i,E2+eta*1i],...
    'ArrayValued',true));

G1 = (2*1i/pi*integral(G1,E1,E2,...
    'Waypoints',[E1+eta*1i,E2+eta*1i],...
    'ArrayValued',true));
G2 = (2*1i/pi*integral(G2,E1,E2,...
    'Waypoints',[E1+eta*1i,E2+eta*1i],...
    'ArrayValued',true));
GE = (2*1i/pi*integral(GE,E1,E2,...
    'Waypoints',[E1+eta*1i,E2+eta*1i],...
    'ArrayValued',true));



% n = arrayfun(dn,Z,'UniformOutput',false);
% n = reshape(cell2mat(n),[Nelem,Npath,Npath]);
% 

% 
% tic;
% Q1 = quadv(dn,E1+eta1*1i,E2+eta1*1i,1e-10);
% toc;
% tic;
% Q2 = quadv(dn,E2+eta1*1i,E2+eta2*1i,1e-10);
% toc;
% tic;
% Q3 = quadv(dn,E2+eta2*1i,E1+eta2*1i,1e-10);
% toc;
% tic;
% Q4 = quadv(dn,E1+eta2*1i,E1+eta1*1i,1e-10);
% toc;
% Qn = [Q1,Q2,Q3,Q4];
% Q = sum(Qn,2);
% 
% 
% n = real(1i*Q/pi);
% 
% plot(x,n)