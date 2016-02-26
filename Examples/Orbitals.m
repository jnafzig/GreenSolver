% Visualize left and right orbitals from shooting method.

clear;
x = [-1/2,1/2]';
Nx = numel(x);

shoot = solver(x);

vL = 0; 
vR = 0;
v = [vL,-2,vR]';

k = @(E,v) sqrt(2*(E-v));

xval = linspace(-10,10,101)';

dens = @(E,v) density(xval,x,k(E,v),shoot(E,v));
densint = @(E,v) densityintegral(x,k(E,v),shoot(E,v));

E0 = min(v)-1;

Nplot = 5;

ER = linspace(-3/2,3/2,Nplot)+1i;
ER(1) = ER(1);

n = zeros(numel(xval),Nplot);
dos = zeros(numel(xval),Nplot);
N = zeros(numel(x)+1,Nplot);

for i = 1:Nplot
    mu = ER(i);
    
    sol = shoot(mu,v);
    
    phiR(:,i) = orbital(xval,x,k(mu,v),sol(:,2));
    phiL(:,i) = orbital(xval,x,k(mu,v),sol(:,1));

    
    dos(:,i) = dens(mu,v);

    
    
%     R = (E0+mu)/2;
%     A = mu-R;
%     E = @(theta) R + A*exp(1i*theta);
%     dEdt = @(theta) 1i*A*exp(1i*theta);
%     
%     tic;
%     n(:,i) = integral(@(theta) dens(E(theta),v)*dEdt(theta),pi,0,...
%                 'ArrayValued',true,...
%                 'RelTol',eps,...
%                 'AbsTol',eps);
%     n(:,i) = (1i*n(:,i)+conj(1i*n(:,i)))/(pi);
%     toc;
% 
% 
%     tic;
%     N(:,i) = integral(@(theta) -densint(E(theta),v)*dEdt(theta),0,pi,...
%                 'ArrayValued',true, ...
%                 'RelTol',eps, ...
%                 'AbsTol',eps);
%     N(:,i) = (1i*N(:,i)+conj(1i*N(:,i)))/(pi);
%     toc;
end


% xstair = [min(xval);x;max(xval)];
% Nstair = [N;N(end,:)];
% vstair = [v;v(end)];
% subplot(2,1,1);
% plot(xval,n); %hold on; stairs(xstair,Nstair,'k--'); hold off;
% xlim([min(xval),max(xval)]);
% 
% title('density');
% 
% subplot(2,1,2);
% plot(xval,dos); 
% 
% xlim([min(xval),max(xval)]);
% 
% title('dos');
plot(xval,[phiL(:,1),phiR(:,1)]); 