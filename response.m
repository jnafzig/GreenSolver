function [ chi] = response( x1, x2, x, k, solution)
    %SOLUTION provide stepwise function for evaluating solutions
    xlt = min(x1,x2);
    xgt = max(x1,x2);
    
    ilt = sum(bsxfun(@gt,xlt,reshape(x,[1,1,numel(x)])),3)+1;
    igt = sum(bsxfun(@gt,xgt,reshape(x,[1,1,numel(x)])),3)+1;
    
    A = solution(1:2:end,1);
    C = solution(1:2:end,2);
    B = solution(2:2:end,1);
    D = solution(2:2:end,2);
    
    W = mean(2.0i*k.*(B.*C-A.*D));
    
    eikxlt = exp(1i*k(ilt).*xlt);
    eikxgt = exp(1i*k(igt).*xgt);
    
    philt = A(ilt).*eikxlt ...
        + B(ilt)./eikxlt;
    phigt = C(igt).*eikxgt ...
        + D(igt)./eikxgt;
    
    chi = philt.^2.*phigt.^2./W.^2;

end
