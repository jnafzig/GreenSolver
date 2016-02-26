% This is an attempt to handle non-zero temperature.  It may be working now

clear;
RelTol = 1e-7;
AbsTol = 1e-7;


x = [-3/2,-1/2,1/2,3/2]';
dx = [0;diff(x);0];
Nx = numel(x);

xval = linspace(-10/2,5/2,101)';
[x1,x2] = ndgrid(xval,xval);    

shoot = solver(x);

vL = 0; 
vR = 0;
v = [vL,-2,0,-2,vR]';

k = @(E,v) sqrt(2*(E-v));

dens = @(E,v) density(xval,x,k(E,v),shoot(E,v));
densint = @(E,v) densityintegral(x,k(E,v),shoot(E,v));

green = @(E,v) greens(x1,x2,x,k(E,v),shoot(E,v));
greenint = @(E,v) greensintegral(x,k(E,v),shoot(E,v));

resp = @(E,v) response(x1,x2,x,k(E,v),shoot(E,v));
respint = @(E,v) responseintegral(x,k(E,v),shoot(E,v));

Emin = min(v)-1;
Emax = 2;

mu = -.1;
beta = 10;

R = (Emin+Emax)/2;
A = Emax-R;
E = @(theta) R + A*exp(1i*theta);
dEdt = @(theta) 1i*A*exp(1i*theta);
f = @(E) 1./(exp(beta*(E-mu))+1);

tic;
n = integral(@(E) dens(E,v)*f(E),Emin,Emax,...
            'WayPoints',[Emin+pi*1i/(2*beta),Emax+pi*1i/(2*beta)],...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
n = (1i*n+conj(1i*n))/(pi);
toc;



tic;
N = integral(@(E) densint(E,v)*f(E),Emin,Emax,...
            'WayPoints',[Emin+pi*1i/(2*beta),Emax+pi*1i/(2*beta)],...
            'ArrayValued',true, ...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
N = (1i*N+conj(1i*N))/(pi);
toc;

tic;
chi = integral(@(E) resp(E,v)*f(E),Emin,Emax,...
            'WayPoints',[Emin+pi*1i/(2*beta),Emax+pi*1i/(2*beta)],...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
chi = 2*(1i*chi+conj(1i*chi))/(pi);
toc;

tic;
chiint = integral(@(E) respint(E,v)*f(E),Emin,Emax,...
            'WayPoints',[Emin+pi*1i/(2*beta),Emax+pi*1i/(2*beta)],...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
chiint = 2*(1i*chiint+conj(1i*chiint))/(pi);
toc;

tic;
densmat = integral(@(E) green(E,v)*f(E),Emin,Emax,...
            'WayPoints',[Emin+pi*1i/(2*beta),Emax+pi*1i/(2*beta)],...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
densmat = (1i*densmat+conj(1i*densmat))/(pi);
toc;

tic;
densmatint = integral(@(E) greenint(E,v)*f(E),Emin,Emax,...
            'WayPoints',[Emin+pi*1i/(2*beta),Emax+pi*1i/(2*beta)],...
            'ArrayValued',true,...
            'RelTol',RelTol,...
            'AbsTol',AbsTol);
densmatint = (1i*densmatint+conj(1i*densmatint))/(pi);
toc;

xstair = [min(xval);x;max(xval)];
Nstair = [N;N(end)];
vstair = [v;v(end)];
subplot(2,2,1);
plot(xval,n);hold on; stairs(xstair,Nstair,'k--'); hold off;
xlim([min(xval),max(xval)]);

title('density and integrated density');

subplot(2,2,3);
stairs(xstair,vstair,'r');

ylim([min(v)-1/2,max(v)+1/2]);
xlim([min(xval),max(xval)]);

title('potential');

subplot(2,2,4);
contourf(x1,x2,chi);
colorbar;

title('response');

subplot(2,2,2);
contourf(x1,x2,densmat); 
colorbar;

title('density matrix');
