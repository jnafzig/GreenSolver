clear;
x = 0;
dx = [0;diff(x);0];
Nx = numel(x);

shoot = solver(x);

vL = -5; 
vR = 0;
v = [vL,vR]';

k = @(E,v) sqrt(2*(E-v));

xval = linspace(-10,5/2,101)';
[x1,x2] = ndgrid(xval,xval);

dens = @(E,v) density(xval,x,k(E,v),shoot(E,v));
densint = @(E,v) densityintegral(x,k(E,v),shoot(E,v));

green = @(E,v) greens(x1,x2,x,k(E,v),shoot(E,v));
greenint = @(E,v) greensintegral(x,k(E,v),shoot(E,v));

resp = @(E,v) response(x1,x2,x,k(E,v),shoot(E,v));
respint = @(E,v) responseintegral(x,k(E,v),shoot(E,v));

E0 = min(v)-1;

Nplot = 5;

muR = linspace(min(v),-.5,Nplot);
muR(1) = muR(1)+.2;

n = zeros(numel(xval),Nplot);
phi = zeros(numel(xval),Nplot);
N = zeros(numel(x)+1,Nplot);

for i = 1:Nplot
    mu = muR(i);
    
    sol = shoot(mu,v);
    phi(:,i) = orbital(xval,x,k(mu,v),sol(:,2))*1e-1 + mu;    
    
    R = (E0+mu)/2;
    A = mu-R;
    E = @(theta) R + A*exp(1i*theta);
    dEdt = @(theta) 1i*A*exp(1i*theta);
    
    tic;
    n(:,i) = integral(@(theta) dens(E(theta),v)*dEdt(theta),pi,0,...
                'ArrayValued',true,...
                'RelTol',eps,...
                'AbsTol',eps);
    n(:,i) = (1i*n(:,i)+conj(1i*n(:,i)))/(pi);
    toc;


    tic;
    N(:,i) = integral(@(theta) -densint(E(theta),v)*dEdt(theta),0,pi,...
                'ArrayValued',true, ...
                'RelTol',eps, ...
                'AbsTol',eps);
    N(:,i) = (1i*N(:,i)+conj(1i*N(:,i)))/(pi);
    toc;
end

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
Nstair = [N;N(end,:)];
vstair = [v;v(end)];
subplot(2,2,1);
plot(xval,n); %hold on; stairs(xstair,Nstair,'k--'); hold off;
xlim([min(xval),max(xval)]);

title('density');

subplot(2,2,3);
stairs(xstair,vstair,'r');hold on; plot(xval,phi); hold off;

ylim([min(v)-1/2,max(v)+1/2]);
xlim([min(xval),max(xval)]);

title('potential + orbitals');

subplot(2,2,4);
contourf(x1,x2,chi);
colorbar;

title('response');

subplot(2,2,2);
contourf(x1,x2,densmat); 
colorbar;

title('density matrix');
