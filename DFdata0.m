function pde = DFdata0
%% DFdata0 data for Darcy-Forchheimer equation 
%
%       - K(|u|) u  + grad p  = f   in \Omega   
%                      div u  = g   in \Omega   
%                   u \dot n  = g_N   on \Gamma   
%
%
%     u = [x+y, x-y];
%     p = x^3 + y^3;
%     f = [(1+beta*sqrt(2x^2 + 2y^2))(x+y)+3x^2, (1+beta*sqrt(2x^2 +
%     2y^2))(x-y)+3y^2]
%     g = 0;
%     g_N = 1+y, x= 1;
%           1-y, x=-1;
%           x-1, y=1;
%           -x-1 y=-1;
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.
beta = 30;
pde = struct('f',@f, 'exactp',@exactp,'exactu',@exactu,'g_D',@g_D,'g_N',@g_N, ...
    'K',@K, 'exactudiv',  @exactudiv, 'beta', beta, 'sigma', @sigma);

    % load data (right hand side function)
    function rhs =  f(p)
        x = p(:,1); y = p(:,2);
        rhs(:, 1) =  (1+beta*sqrt(2*x.^2 + 2*y.^2)) .*(x+y) + 3*x.^2;
        rhs(:, 2) =  (1+beta*sqrt(2*x.^2 + 2*y.^2)) .*(x-y) + 3*y.^2;

%         rhs(:, 1) =  0*x + 2;
%         rhs(:, 2) =  0*y + 2; 
    end
    function s = K(unorm)
       s = beta*(unorm) + 1;
    end

    function s = sigma(unorm)
       s = beta*(unorm);
    end
    % exact solution
    function s =  exactp(pt)
        x = pt(:,1); y = pt(:,2);
        s = x.^3 + y.^3;
%         s = 0.*x;
    end
    % the derivative of the exact solution
    function s = exactu(p)
        x = p(:,1); y = p(:,2);
        s(:,1) = x+y;
        s(:,2) = x-y;
% 
%         s(:,1) = 0.*x + 2;
%         s(:,2) = 0.*y + 2;
    end
    function u = g_D(p)
        u = exactudiv(p);
    end
    function u = exactudiv(p)
        x = p(:,1); y = p(:,2);
        u = 0.*x;
    end
    % Neumann boundary condition 
    function f = g_N(p,vargin)
        if nargin > 1
            f = dot(exactu(p),vargin,2);
        else
            f = zeros(size(p,1),1);
            x = p(:,1); y = p(:,2);
            uprime = exactu(p);
            leftbd = (abs(x+1)<eps);  % n = (-1,0); 
            f(leftbd) = - uprime(leftbd,1);
            rightbd = (abs(x-1)<eps); % n = (1,0); 
            f(rightbd) = uprime(rightbd,1);
            topbd = (abs(y-1)<eps);   % n = (0,1)
            f(topbd) = uprime(topbd,2);
            bottombd = (abs(y+1)<eps);% n = (0,-1)
            f(bottombd) = - uprime(bottombd,2);    
        end
    end
end