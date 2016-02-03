% Test Script for playing around with finite difference matrices

clear;

Nelem = 501;
Nstencil = 7;
locR = 0:Nstencil;
locL = -Nstencil:0;

coeff = fd_coeff(locR,2,1);
D2R = spdiags(repmat(coeff',[Nelem,1]), locR,Nelem,Nelem);
bD2R = spdiags(repmat(coeff(2:end)',[Nstencil,1]), locR(2:end)-Nstencil,Nstencil,Nstencil);

coeff = fd_coeff(locR,1,1);
D1R = spdiags(repmat(coeff',[Nelem,1]), locR,Nelem,Nelem);
bD1R = spdiags(repmat(coeff(2:end)',[Nstencil,1]), locR(2:end)-Nstencil,Nstencil,Nstencil);

coeff = fd_coeff(locL,2,1);
D2L = spdiags(repmat(coeff',[Nelem,1]), locL,Nelem,Nelem);
bD2L = spdiags(repmat(coeff(1:end-1)',[Nstencil,1]), locL(1:end-1)+Nstencil,Nstencil,Nstencil);

coeff = fd_coeff(locL,1,1);
D1L = spdiags(repmat(coeff',[Nelem,1]), locL,Nelem,Nelem);
bD1L = spdiags(repmat(coeff(1:end-1)',[Nstencil,1]), locL(1:end-1)+Nstencil,Nstencil,Nstencil);


xfull = linspace(-10,10,Nelem)';
x = xfull(1:Nelem);
dx = x(2)-x(1);

R = -1;
L = 4;
v0 = -3;

v = 0*x;
v(x>=(R-L/2)&x<=(R+L/2)) = v0;

HR = sparse(-1/2*D2R/dx^2+diag(v));
bHR = -1/2*bD2R/dx^2;
HL = sparse(-1/2*D2L/dx^2+diag(v));
bHL = -1/2*bD2L/dx^2;

E = -1;
k = sqrt(2*E);
    
bx = (1:Nstencil)'*dx;
bphi = exp(1i*k*bx);
bHval = [zeros(Nelem-Nstencil,1);bHR*bphi];
        

phi = (HR-E*speye(Nelem))\(-bHval);
    



% [phiR,dphiR] = shoot(E,HR,bHR,D1R,bD1R,dx, true);
% [phiL,dphiL] = shoot(E,HL,bHL,D1L,bD1L,dx, false);
% G = phiR.*phiL/mean(phiL.*dphiR-phiR.*dphiL);
% 
% W = (phiL.*dphiR-phiR.*dphiL);
% plot(x,[real(G),imag(G)])

