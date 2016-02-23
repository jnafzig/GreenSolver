function [ N] = densityintegral( x, k, solution)
        
    Nx = numel(x);
    i = reshape([1:Nx;(1:Nx)+1],[2*Nx,1]);
    j = reshape([1:Nx;(1:Nx)],[2*Nx,1]);
    
    A = solution(1:2:end,1);
    B = solution(2:2:end,1);
    C = solution(1:2:end,2);
    D = solution(2:2:end,2);
    
    
    AC = A.*C;
    BCAD = B.*C + A.*D;
    BD = B.*D;
    
    W = mean(2.0i*k.*(B.*C-A.*D));
    
    N = (BCAD(i).*x(j) ...
        + BD(i).*exp(-2.0i*k(i).*x(j))./(-2.0i*k(i)) ...
        + AC(i).*exp(2.0i*k(i).*x(j))./(2.0i*k(i)) ) / W;
        
    N = [N(1:2:end);0]-[0;N(2:2:end)];
end
  
 