function testresponse( testCase )
    %TESTRESPONSE test whether response correctly predicts density
    % response to small change in potential 
    
    x = [-1,0,1]';

    shoot = solver(x);

    
    % testing set of potentials
    vtest = [[ 0     7     1    0]; ...
             [-3     8     5    0]; ...
             [ 0     2     2   -4]; ...
             [10     8    -9   -8]; ...
             [-4     9    -1   -1]];
    Ntest = size(vtest,1);
         
    dv = sqrt(eps);
    vdiff = dv*[2  0  3  -3]'; % small change in the potential
    
    for i = 1:Ntest
        
        v1 = vtest(i,:)';
        v2 = v1+vdiff;

        k = @(E,v) sqrt(2*(E-v));

        densint = @(E,v) densityintegral(x,k(E,v),shoot(E,v));
        respint = @(E,v) responseintegral(x,k(E,v),shoot(E,v));

        E0 = min(v1)-1;
        mu = -.5;

        R = (E0+mu)/2;
        A = mu-R;
        E = @(theta) R + A*exp(1i*theta);
        dEdt = @(theta) 1i*A*exp(1i*theta);

        N1 = integral(@(theta) -densint(E(theta),v1)*dEdt(theta),0,pi,...
                    'ArrayValued',true, ...
                    'RelTol',eps, ...
                    'AbsTol',eps);
        N1 = (1i*N1+conj(1i*N1))/(pi);

        N2 = integral(@(theta) -densint(E(theta),v2)*dEdt(theta),0,pi,...
                    'ArrayValued',true,...
                    'RelTol',eps,...
                    'AbsTol',eps);
        N2 = (1i*N2+conj(1i*N2))/(pi);

        chi = integral(@(theta) -respint(E(theta),v1)*dEdt(theta),0,pi,...
                    'ArrayValued',true,...
                    'RelTol',eps,...
                    'AbsTol',eps);
        chi = 2*(1i*chi+conj(1i*chi))/(pi);

        % The change in density should be equal to response * change in
        % potential
        Check = (N2-N1) - chi*vdiff;

        testCase.verifyEqual(max(abs(Check)),0,'AbsTol',1e2*eps);
    end
end

