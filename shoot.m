function [ phi,dphi ] = shoot( E, H, bH, D1, bD1, bV, dx, right )
    % SHOOT( E, H, bH, x, bx )
    
    [Nelem,~] = size(H);
    [bNelem,~] = size(bH);
    k = sqrt(2*(E-bV));
    
    if right
        bx = (1:bNelem)'*dx;
        if imag(k) > 0
            bphi = exp(1i*k*bx);
        else
            bphi = exp(-1i*k*bx);
        end
        bHval = [zeros(Nelem-bNelem,1);bH*bphi];
        bD1val = [zeros(Nelem-bNelem,1);bD1*bphi];
    else % left
        bx = (-bNelem:-1)'*dx;
        if imag(k) > 0
            bphi = exp(-1i*k*bx);
        else
            bphi = exp(1i*k*bx);
        end
        bHval = [bH*bphi;zeros(Nelem-bNelem,1)];
        bD1val = [bD1*bphi;zeros(Nelem-bNelem,1)];
    end
    
    phi = (H-E*speye(Nelem))\(-bHval);
    
    dphi = (D1*phi + bD1val)/dx;

end

