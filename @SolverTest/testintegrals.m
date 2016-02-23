function testintegrals( testCase )
    %TEST test that the density integral function correctly integrates the
    % density
    x = 0;
    
    kval = [ 1i; 1i ];

    sol = [ ...
      [1 1];...
      [2 1];...
      [0 2];...
      [1 0] ];
 
    Check1 = densityintegral(x,kval,sol);
    Check2 = greensintegral(x,kval,sol);
    Check3 = responseintegral(x,kval,sol);

    testCase.verifyEqual(Check1,[-1/6; 0],'AbsTol',1e2*eps);
    testCase.verifyEqual(Check2,[[0 -2/3];[-2/3 2/3]],'AbsTol',1e2*eps);
    testCase.verifyEqual(Check3,[[1/4 1/3];[1/3 -1/9]],'AbsTol',1e2*eps);
     
end

