clear;

Nelem = 101;
Nstencil = 7;

x = linspace(-10,10,Nelem)';
dx = x(2)-x(1);

X0 = 0;
vL = -3;
vR = 0;

v = 0*x;
v(x<X0) = vL;
v(x>=X0) = vR;

[G,Wv] = green(Nelem,Nstencil,dx);
dn = @(E) G(E,v);
W = @(E) Wv(E,v);

E1 = -4;
E2 = 1;
eta1 = 1;
eta2 = -1;

Npath = 100;

[ZR,ZI] = ndgrid(linspace(E1,E2,Npath),linspace(eta2,eta1,Npath));
Z = ZR+ZI*1i;

integrand = arrayfun(dn,Z,'UniformOutput',false);
integrand = reshape(cell2mat(integrand),[Nelem,Npath,Npath]);

tic;
Q1 = quadv(dn,E1+eta1*1i,E2+eta1*1i,1e-10);
toc;
tic;
Q2 = quadv(dn,E2+eta1*1i,E2+eta2*1i,1e-10);
toc;
tic;
Q3 = quadv(dn,E2+eta2*1i,E1+eta2*1i,1e-10);
toc;
tic;
Q4 = quadv(dn,E1+eta2*1i,E1+eta1*1i,1e-10);
toc;


Qn = [Q1,Q2,Q3,Q4];
Q = sum(Qn,2);


n = real(1i*Q/pi);

plot(x,n)

Wv = arrayfun(W,Z);

