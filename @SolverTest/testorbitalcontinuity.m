function testorbitalcontinuity( testCase )
    %TEST test orbital function and continuity simultaneously
    
    x = [-2;-1;0;1;2];
    
    shoot = solver(x);
    
    % testing set of potentials
    vtest = [[ 0     7     1     3     7     0]; ...
             [-3     8     5     1    -2    0]; ...
             [ 0     2     2    -2     2     -4]; ...
             [10     8    -9     9     8    -8]; ...
             [-4     9    -1     7     9    -1]];
         
    Ntest = size(vtest,1);
         
    for i = 1:Ntest
        v = vtest(i,:)';
        E = -1/2+1i;
        k = sqrt(2*(E-v));
        sol = shoot(E,v);              % find left and right orbitals

        Check1 = orbital([x,x+eps],x,k,sol); % evaluate orbital at and just to
                                             % right of grid points

        Check2 = orbital([x,x+eps],x,k,sol,1); % evaluate derivative of orbital 
                                         % at and just to right of grid points

        testCase.verifyEqual(max(abs(diff(Check1,1,2))),0,'AbsTol',1e-10);
        testCase.verifyEqual(max(abs(diff(Check2,1,2))),0,'AbsTol',1e-10);
    end
end

