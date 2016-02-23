function [ G] = greens( x1, x2, x, k, solution)
    %SOLUTION provide stepwise function for evaluating solutions
    xlt = min(x1,x2);
    xgt = max(x1,x2);
    
    i1 = sum(bsxfun(@gt,xlt,reshape(x,[1,1,numel(x)])),3)+1;
    i2 = sum(bsxfun(@gt,xgt,reshape(x,[1,1,numel(x)])),3)+1;    
   
    A = solution(1:2:end,1);
    B = solution(2:2:end,1);
    C = solution(1:2:end,2);
    D = solution(2:2:end,2);
    
    W = mean(2.0i*k.*(B.*C-A.*D));
    
    eikx1 = exp(1i*k(i1).*xlt);
    eikx2 = exp(1i*k(i2).*xgt);
    G = (A(i1).*eikx1+B(i1)./eikx1)...
      .*(C(i2).*eikx2+D(i2)./eikx2)/W;

    
end
