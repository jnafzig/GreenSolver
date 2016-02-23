% Example script for a simple constant finite well

clear;
x = [-1/2,1/2]';
dx = [0;diff(x);0];
Nx = numel(x);

shoot = solver(x);

vL = 0; 
vR = 0;
v = [vL,-20,vR]';

k = @(E,v) sqrt(2*(E-v));

xval = linspace(-3/2,3/2,101)';
[x1,x2] = ndgrid(xval,xval);
    
dens = @(E,v) density(xval,x,k(E,v),shoot(E,v));
densint = @(E,v) densityintegral(x,k(E,v),shoot(E,v));

green = @(E,v) greens(x1,x2,x,k(E,v),shoot(E,v));
greenint = @(E,v) greensintegral(x,k(E,v),shoot(E,v));

resp = @(E,v) response(x1,x2,x,k(E,v),shoot(E,v));
respint = @(E,v) responseintegral(x,k(E,v),shoot(E,v));

E0 = min(v)-1;
mu = -.5;

R = (E0+mu)/2;
A = mu-R;
E = @(theta) R + A*exp(1i*theta);
dEdt = @(theta) 1i*A*exp(1i*theta);

tic;
n = integral(@(theta) dens(E(theta),v)*dEdt(theta),pi,0,...
            'ArrayValued',true,...
            'RelTol',eps,...
            'AbsTol',eps);
n = (1i*n+conj(1i*n))/(pi);
toc;


tic;
N = integral(@(theta) -densint(E(theta),v)*dEdt(theta),0,pi,...
            'ArrayValued',true, ...
            'RelTol',eps, ...
            'AbsTol',eps);
N = (1i*N+conj(1i*N))/(pi);
toc;

tic;
chi = integral(@(theta) resp(E(theta),v)*dEdt(theta),pi,0,...
            'ArrayValued',true,...
            'RelTol',eps,...
            'AbsTol',eps);
chi = 2*(1i*chi+conj(1i*chi))/(pi);
toc;

tic;
chiint = integral(@(theta) -respint(E(theta),v)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',eps,...
            'AbsTol',eps);
chiint = 2*(1i*chiint+conj(1i*chiint))/(pi);
toc;

tic;
densmat = integral(@(theta) green(E(theta),v)*dEdt(theta),pi,0,...
            'ArrayValued',true,...
            'RelTol',eps,...
            'AbsTol',eps);
densmat = (1i*densmat+conj(1i*densmat))/(pi);
toc;

tic;
densmatint = integral(@(theta) -greenint(E(theta),v)*dEdt(theta),0,pi,...
            'ArrayValued',true,...
            'RelTol',eps,...
            'AbsTol',eps);
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
