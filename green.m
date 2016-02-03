function [ GE, WE, GExx ] = green( Nelem,Nstencil,dx )
    
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

    HR = -1/2*D2R/dx^2;
    bHR = -1/2*bD2R/dx^2;
    HL = -1/2*D2L/dx^2;
    bHL = -1/2*bD2L/dx^2;
    
    GE = @(E,v) G(E,v);
    WE = @(E,v) W(E,v);
    GExx = @(E,v) Gxx(E,v);
    
    function G = G(E,v)

        bVR = v(end);
        bVL = v(1);
        
        [phiR,dphiR] = shoot(E,HR+spdiags(v,0,Nelem,Nelem),bHR,D1R,bD1R,bVR,dx, true);
        [phiL,dphiL] = shoot(E,HL+spdiags(v,0,Nelem,Nelem),bHL,D1L,bD1L,bVL,dx, false);
        
        W = (phiL.*dphiR-phiR.*dphiL);
        G = phiR.*phiL./mean(W);

    end
    
    function W = W(E,v)

        bVR = v(end);
        bVL = v(1);
        
        [phiR,dphiR] = shoot(E,HR+spdiags(v,0,Nelem,Nelem),bHR,D1R,bD1R,bVR,dx, true);
        [phiL,dphiL] = shoot(E,HL+spdiags(v,0,Nelem,Nelem),bHL,D1L,bD1L,bVL,dx, false);
        
        W = mean(phiL.*dphiR-phiR.*dphiL);

        
    end
    
    function G = Gxx(E,v)

        bVR = v(end);
        bVL = v(1);
        
        [phiR,dphiR] = shoot(E,HR+spdiags(v,0,Nelem,Nelem),bHR,D1R,bD1R,bVR,dx, true);
        [phiL,dphiL] = shoot(E,HL+spdiags(v,0,Nelem,Nelem),bHL,D1L,bD1L,bVL,dx, false);
        
        W = mean(phiL.*dphiR-phiR.*dphiL);

        G = (triu(phiL*transpose(phiR),1)+tril(phiR*transpose(phiL)))...
            / mean(W);
        
    end
    
end

