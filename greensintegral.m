function [ G] = greensintegral( x, k, solution)
    %SOLUTION provide stepwise function for evaluating solutions
    Nx = numel(x);
    i = reshape([1:Nx;(1:Nx)+1],[2*Nx,1]);
    j = reshape([1:Nx;(1:Nx)],[2*Nx,1]);
    
    [h1,h2] = ndgrid(1:2*Nx,1:2*Nx);
    
    hlt = min(h1,h2);
    hgt = max(h1,h2);

    A = solution(1:2:end,1);
    B = solution(2:2:end,1);
    C = solution(1:2:end,2);
    D = solution(2:2:end,2);

    W = mean(2.0i*k.*(B.*C-A.*D));
    
    Ilt = -1i*(A(i).*exp(1i*k(i).*x(j)) ...
            - B(i).*exp(-1i*k(i).*x(j))) ./k(i);
    Igt = -1i*(C(i).*exp(1i*k(i).*x(j)) ...
            - D(i).*exp(-1i*k(i).*x(j))) ./k(i);
    
    G = Ilt(hlt).*Igt(hgt)/W ;
    G = [G(1:2:end,:);zeros(1,2*Nx)]-[zeros(1,2*Nx);G(2:2:end,:)];
    G = [G(:,1:2:end),zeros(Nx+1,1)]-[zeros(Nx+1,1),G(:,2:2:end)];

    Gdiag = x(j)./k(i).^2;
    Gdiag = [Gdiag(1:2:end);0]-[0;Gdiag(2:2:end)];

    G = G + diag(Gdiag);
end
