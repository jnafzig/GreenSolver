function [ solver_fh ] = solver( x )
    %SOLVER 
    
    Nx = numel(x);

    r = @(N) (1:2:2*N)';

    Mi = [[r(Nx);r(Nx);r(Nx)+1;r(Nx)+1];[r(Nx);r(Nx);r(Nx)+1;r(Nx)+1]];
    Mj = [[r(Nx);r(Nx)+1;r(Nx);r(Nx)+1]+2;[r(Nx);r(Nx)+1;r(Nx);r(Nx)+1]];
    Mval = @(k) [exp(k(2:end).*x*1i);...
                exp(-k(2:end).*x*1i);...
                1i*k(2:end).*exp(k(2:end).*x*1i);...
                -1i*k(2:end).*exp(-k(2:end).*x*1i);...
                -exp(k(1:end-1).*x*1i);...
                -exp(-k(1:end-1).*x*1i);... 
                -1i*k(1:end-1).*exp(k(1:end-1).*x*1i);...
                1i*k(1:end-1).*exp(-k(1:end-1).*x*1i)];

            
    M = @(k) sparse(Mi,Mj,Mval(k),2*Nx,2*(Nx+1));
    bcL = [[0;1];zeros(2*Nx,1)];
    bcR = [zeros(2*Nx,1);[1;0]];
    bcMatchL = [eye(2),zeros(2,2*Nx)];
    bcMatchR = [zeros(2,2*Nx),eye(2)];
    
    solver_fh = @(E,v) shoot(E,v);
    
    function solution = shoot(E,v)
        k = sqrt(2*(E-v));
        Match = M(k);
        solution = [[bcMatchL;Match]\bcL,...
                   [Match;bcMatchR]\bcR];
               
        solution(1:2,1) = [0;1];
        solution(end-1:end,2) = [1;0];

    end
end

