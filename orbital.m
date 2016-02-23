function [ phi] = orbital( xval, x, k, solution, order)
    %SOLUTION provide stepwise function for evaluating solutions
    if nargin < 5
        order = 0;
    end
    
    i = sum(bsxfun(@gt,xval,reshape(x,[1,1,numel(x)])),3)+1;
    
    A = solution(1:2:end);
    B = solution(2:2:end); 
    
    phi = (1i*k(i)).^order.*A(i).*exp(1i*k(i).*xval) + ...
           (-1i*k(i)).^order.*B(i).*exp(-1i*k(i).*xval);

end