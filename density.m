function [ dn] = density( xval, x, k, solution)
    %SOLUTION provide stepwise function for evaluating solutions
    [Nx,Mx] = size(xval);
    xval = reshape(xval,[numel(xval),1]);
    
    i = sum(bsxfun(@gt,xval,reshape(x,[1,1,numel(x)])),3)+1;
    
    A = solution(1:2:end,1);
    B = solution(2:2:end,1);
    C = solution(1:2:end,2);
    D = solution(2:2:end,2);
    
    AC = A.*C;
    BCAD = B.*C + A.*D;
    BD = B.*D;
    
    W = mean(2.0i*k.*(B.*C-A.*D));
    
    dn = (BCAD(i) ...
        + BD(i).*exp(-2.0i*k(i).*xval) ...
        + AC(i).*exp(2.0i*k(i).*xval))/W;
    
    dn = reshape(dn,[Nx,Mx]);
end

 