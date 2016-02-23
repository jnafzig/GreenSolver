function [ chi] = responseintegral( x, k, solution)
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
   
    Igt2 = C(i).^2.*exp(2i*k(i).*x(j))./(2i*k(i)) ...
        + D(i).^2.*exp(-2i*k(i).*x(j))./(-2i*k(i))...
        + 2*C(i).*D(i).*x(j);
    
    Ilt2 = A(i).^2.*exp(2i*k(i).*x(j))./(2i*k(i)) ...
        + B(i).^2.*exp(-2i*k(i).*x(j))./(-2i*k(i)) ...
        + 2*A(i).*B(i).*x(j);
    
    chi = Ilt2(hlt).*Igt2(hgt)/W^2;
    
    chidiag = ( C(i).*A(i).*(1 - 1i*k(i).*x(j)).*exp(2i*k(i).*x(j)) ...
              - D(i).*B(i).*(1 + 1i*k(i).*x(j))./exp(2i*k(i).*x(j)) ...
              + 1i*(A(i).*D(i) + C(i).*B(i)).*k(i).*x(j) ) ...
              ./(4*(A(i).*D(i) - C(i).*B(i)).*k(i).^4);
  
    chidiag = [chidiag(1:2:end);0]-[0;chidiag(2:2:end)];
    
    chi = [chi(1:2:end,:);zeros(1,2*Nx)]-[zeros(1,2*Nx);chi(2:2:end,:)];
    chi = [chi(:,1:2:end),zeros(Nx+1,1)]-[zeros(Nx+1,1),chi(:,2:2:end)];
    
    chi = chi + diag(chidiag);
end
