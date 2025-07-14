function pde = DFdata1
%% DFdata0 data for Darcy-Forchheimer equation 
%
%       - K(|u|) u  + grad p  = f   in \Omega   
%                      div u  = g   in \Omega   
%                   u \dot n  = g_N   on \Gamma   
%
%
%     u = [(x+1)^2/4, (x+1)(y+1)/2];
%     p = x^3 + y^3;
%     f = [(x+1)^2/4*(1+beta*(x+1)/4*sqrt((x+1)^2 + 4(y+1)^2))+3x^2, 
%          -(x+1)(y+1)/2*(1+beta*(x+1)/4*sqrt((x+1)^2 + 4(y+1)^2))+3y^2]
%     g = 0;
%     Dp = (2*pi*cos(2*pi*x)*cos(2*pi*y), -2*pi*sin(2*pi*x)*sin(2*pi*y));
%     g_N = 1, x= 1;
%           0, x=-1;
%           -x-1, y=1;
%           0, y=-1;
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.
beta = 30;
pde = struct('f',@f, 'exactp',@exactp,'exactu',@exactu,'g_D',@g_D,'g_N',@g_N, ...
    'K',@K, 'exactudiv',  @exactudiv, 'Ju', @Ju, 'Juinv', @Juinv, 'beta', beta,'sigma', @sigma);

    % load data (right hand side function)
    function rhs =  f(p)
        x = p(:,1); y = p(:,2);
        rhs(:, 1) = (x+1).^2/4.*(1+beta*(x+1)/4.*sqrt((x+1).^2 + 4*(y+1).^2)) + 3*x.^2;
        rhs(:, 2) = -(x+1).*(y+1)/2.*(1+beta*(x+1)/4.*sqrt((x+1).^2 + 4*(y+1).^2)) + 3*y.^2;

%         rhs(:, 1) =  0*x + 2;
%         rhs(:, 2) =  0*y + 2; 
    end
    function s = K(unorm)
       s = beta*(unorm) + 1;
    end

    function s = sigma(unorm)
       s = beta*(unorm) + 1;
    end

    function M = Ju(u)
        Nu = length(u);
        NT = Nu/2;
        unorm = sqrt(u(1:NT).^2 + u(NT+1:end).^2);
        s = beta*(unorm) + 1;
        ii = [1:Nu, 1:Nu, 1:NT, NT+1:Nu];
        jj = [1:Nu, 1:Nu, NT+1:Nu, NT+1:Nu];
        ss = [s; s; beta*u.^2./[unorm; unorm]; beta*u(1:NT).*u(NT+1:end)./unorm; ...
            beta*u(1:NT).*u(NT+1:end)./unorm]
        M = sparse(ii, jj, ss, Nu, Nu);
%         M = sparse(1:Nu, 1:Nu, [s; s], Nu, Nu);
%         M = M + sparse(1:Nu, 1:Nu, beta*u.^2./[unorm; unorm], Nu, Nu);
%         M = M + sparse(1:NT, NT+1:Nu,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
%         M = M + sparse( NT+1:Nu, 1:NT,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
    end

    function M = Juinv(u)
        Nu = length(u);
        NT = Nu/2;
        unorm = sqrt(u(1:NT).^2 + u(NT+1:end).^2);
        s = beta*(unorm) + 1;
        detu = (beta*(u(1:NT).^2)./unorm+s).*(beta*(u(NT+1:end).^2)./unorm+s) ...
            - (beta^2)*(u(1:NT).^2).*(u(NT+1:end).^2)./(unorm.^2);
        detuinv = 1./detu;
        M = sparse(1:Nu, 1:Nu, [s + beta*(u(NT+1:end).^2./unorm) ; beta*(u(1:NT).^2)./unorm+s], Nu, Nu);
        M = M + sparse(1:NT, NT+1:Nu,-beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
        M = M + sparse( NT+1:Nu, 1:NT,-beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
        M = [ detuinv; detuinv].*M;
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
        s(:,1) = (x+1).^2/4;
        s(:,2) = -(x+1).*(y+1)/2;
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