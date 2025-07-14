function [p,u,eqn,info] = DFP0P1TPD(node,elem,bdFlag,pde,option)

if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

N = size(node,1); 

maxIt = 2e2;
tol = 1e-9;
flag = 0;

%% Data structure
[Dlambda,area,~] = gradbasis(node,elem);
[elem2edge,edge, ~] = dofedge(elem);
isNeumann(elem2edge(bdFlag(:)==2)) = true;
Neumann = edge(isNeumann,:); 

NT = size(elem,1);
Nu = 2*NT; Np = N; Ndof = Nu + Np;

%% fixed Diffusion coefficient
if ~isfield(pde,'K')
    pde.K = []; 
    K = [];
end
if ~isempty(pde.K) && isnumeric(pde.K)
    K = pde.K;                                 % d is an array
end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end

%% Assemble matrix for divergence operator
d1 = Dlambda(:,:,1).*[area, area];
d2 = Dlambda(:,:,2).*[area, area];
d3 = Dlambda(:,:,3).*[area, area];
Dx = sparse(elem(:), repmat((1:NT)',3,1),...
            [d1(:,1); d2(:,1); d3(:,1)], Np, NT);
Dy = sparse(elem(:), repmat((1:NT)',3,1),...
            [d1(:,2); d2(:,2); d3(:,2)], Np, NT);
B = [Dx Dy];
clear d1 d2 d3

% Assemble right hand side.
fu = zeros(Nu,1);% the right hand side of u
if ~isfield(pde,'f') || (isfield(pde,'f') && isreal(pde.f) && all(pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order is 3
end
if isfield(pde,'f') && ~isempty(pde.f)
    [lambda,w] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
    fp = zeros(NT, 2);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ... 
            + lambda(p,3)*node(elem(:,3),:);
        fp = fp + w(p)*pde.f(pxy);
    end
    fp = fp.*repmat(area,1,2);
    fu = fp(:);
end
normf = norm(fu);
F(1:Nu, 1) = fu;

g = zeros(Np,1); % the right hand side of p
% if ~isfield(pde,'g') || (isreal(pde.g) && (pde.g==0))
%     pde.g = [];
% end
% if ~isfield(option,'gquadorder')
%     option.gquadorder = 2;   % default order
% end
% if ~isempty(pde.g)
% 	[lambda,weight] = quadpts(option.gquadorder);
% 	nQuad = size(lambda,1);
%     for k = 1:nQuad
% 		% quadrature points in the x-y coordinate
% 		pxy = lambda(k,1)*node(elem(:,1),:) ...
% 			+ lambda(k,2)*node(elem(:,2),:) ...
% 			+ lambda(k,3)*node(elem(:,3),:);
% 		gp = pde.g_D(pxy);
% 		g = g + gp*weight(k); % div u = g;
%     end
%     g = g.*area;
% end
% clear gp


gN = zeros(Np,1); % the right hand sife from integration by part of divergence.
if ~isfield(pde,'g_N') || (isreal(pde.g_N) && (pde.g_N==0))
    pde.g_N = [];
end
if ~isfield(option,'gNquadorder')
    option.gNquadorder = 2;   % default order
end
if ~isempty(pde.g_N)
    % modify the rhs to include boundary condition 
    ve = node(Neumann(:,1),:)-node(Neumann(:,2),:);
    edgeLength = sqrt(sum(ve.^2,2)); 
    % compute the integral int_e g_N ds
    [lambda,weight] = quadpts1(option.gNquadorder);
    nQuad = size(lambda,1);
    gNp = zeros(size(Neumann));
    for ip = 1:nQuad
    	pxy = lambda(ip,1)*node(Neumann(:,1),:)+...
              lambda(ip,2)*node(Neumann(:,2),:);               
        gNp(:, 1) = gNp(:, 1) + weight(ip)*pde.g_N(pxy).*edgeLength*lambda(ip,1);
        gNp(:, 2) = gNp(:, 2) + weight(ip)*pde.g_N(pxy).*edgeLength*lambda(ip,2);
    end
    gN = accumarray(Neumann(:, 1), gNp(:, 1), [Np, 1]);
    gN = gN + accumarray(Neumann(:, 2), gNp(:, 2), [Np, 1]);
end
F((Nu+1):Ndof,1) = g + gN;

normg = norm(F((Nu+1):Ndof,1));

%% Initial values
freep = 1:Np-1;
% 
% uoldAll = zeros(Nu,1);
% poldAll = zeros(Np,1);

uoldAll = rand(Nu,1);
poldAll = rand(Np,1);
poldAll(end) = 0;

% Assemble Mass Matrix
if ~isempty(pde.K) && ~isnumeric(pde.K)       % d is a function   
    unorm = sqrt(uoldAll(1:NT).^2 + uoldAll(NT+1:end).^2);
    K = pde.K(unorm); % K is a scalar
end

% M. Mass matrix for P0 element
% M = sparse(1:Nu, 1:Nu, [K.*area;K.*area], Nu, Nu);

M = [K.*area;K.*area];
Minv = 1./M;


% Jinv = pde.Juinv(uoldAll);
% MJinv = [1./area; 1./area].*Jinv;
% S0 = eqn.Bff*(MJinv*(eqn.Bfft));

eqn = struct('M',M,'B',B,...
             'f',F(1:Nu),'g',F(Nu+1:end),...
             'Bff',B(freep,:), 'Bfft',B(freep,:)', 'gf',F(Nu+freep));
Bff  = B(freep,:);
Bfft = Bff';
% S0 = eqn.Bff*(M\(eqn.Bfft));
S0 = eqn.Bff*(Minv.*(eqn.Bfft));

% set up multilevel structure for mg solver
setupOption.solver = 'NO';
setupOption.freeDof =  [true(Np-1, 1); false];
bigS0 = [S0 sparse(Np-1,1); sparse(1,Np-1) 1];
[~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(bigS0,poldAll,elem,setupOption);

% A = [M (B(1:Np-1, :))'; B(1:Np-1, :) sparse(Np-1, Np-1)]; % direct solver
% bigu = A\F(1:end-1);
% unewAll = bigu(1:Nu);
% pnewAll = bigu(Nu+1:end);
% pnewAll(Np) = 0;


%% Error
resu = zeros(maxIt,1);
resp = zeros(maxIt,1);

%% Solve the saddle point problem
stepsize = 0.4;
scaleu = 0.8;
gamma0 = 2;
Vcycle_num = 0;
Bp =  Bfft*poldAll(freep);
tic
for ite = 1:maxIt
%      fprintf('Iteration: %i \n', ite)
    
%     Jinv = pde.Juinv(uoldAll);
%     MJinv = [1./area; 1./area].*Jinv;
%     

     % explicit TPD   
     du_tmp = M.*uoldAll + Bp - eqn.f;
%      du = du_tmp; 
     du = du_tmp./M;
%      du = MJinv*du_tmp;

     unew_tmp = uoldAll - du/scaleu;

     dp_tmp = eqn.gf - eqn.Bff*unew_tmp;
%      dp = S0\dp_tmp;
%      pnew = poldAll(freep) - stepsize*dp;

%     Minv = 1./ spdiags(M, 0);
%     S0 = eqn.Bff*(Minv.*(eqn.Bfft));
%     S = eqn.Bff*(MJinv*(eqn.Bfft));


  % S0 = (1-(stepsize*gamma0))*S0 + (stepsize*gamma0)*S; % explicit IQ
% 
      if mod(ite, 5) == 0
            Minv = 1./M;
            S = eqn.Bff*(Minv.*(eqn.Bfft));
%           S = eqn.Bff*(M\(eqn.Bfft));
          S0 = (S0 + (stepsize*gamma0)*S)/(1+(stepsize*gamma0)); % implicit IQ
          bigS0 = [S0 sparse(Np-1,1); sparse(1,Np-1) 1];
          [Ai, Bi, BBi] = mgmat(S0, Ai, Res, Pro);
      end
%       S = eqn.Bff*(M\(eqn.Bfft));
%       S0 = (S0 + (stepsize*gamma0)*S)/(1+(stepsize*gamma0)); % implicit IQ
%      bigS0 = [S0 zeros(Np-1,1);zeros(1,Np-1) 1];
%      [Ai, Bi, BBi] = mgmat(S0, Ai, Res, Pro);

     [dpAll, info] = mg(bigS0, [dp_tmp;0], elem, option.mg,Ai,Bi,BBi,Res,Pro,isFreeDof);
%      [dpAll, info] = mg(bigS0, [dp_tmp;0], elem, option.mg);
     Vcycle_num  = Vcycle_num  + info.itStep;
     pnew = poldAll(freep) - (stepsize/scaleu)*dpAll(freep);
  
     pnewAll =  poldAll;   % boundary values
     pnewAll(freep) = pnew;

     unewAll = (1-stepsize)*uoldAll + stepsize*unew_tmp;
     
       % M. Mass matrix for P0 element

    if ~isempty(pde.K) && ~isnumeric(pde.K)       % d is a function   
        unorm = sqrt(unewAll(1:NT).^2 + unewAll(NT+1:end).^2);
        K = pde.K(unorm); % K is a scalar
    end
%     M = sparse(1:Nu, 1:Nu, [K.*area;K.*area], Nu, Nu);
    M = [K.*area;K.*area];
    eqn.M = M;

  
    % compute error
    Bp =  Bfft*pnew;
    resu(ite) = norm(M.*unewAll + Bp- eqn.f);
    resp(ite) = norm(eqn.Bff*unewAll  - eqn.gf);
     
    if ite> 1 && max(resu(ite)/normf, resp(ite)/normg) < tol
        flag = 1;
        break;
    end

    uoldAll = unewAll;
    poldAll = pnewAll;
end
time = toc;

% compute L2 norm of p, normalize
if ~isfield(option,'pquadorder')
    option.pquadorder = 2;   % default order
end
[lambda,weight] = quadpts(option.pquadorder);
nQuad = size(lambda,1);
pbar = zeros(NT,1);
for p = 1:nQuad
     pt = pnewAll(elem(:,1))*lambda(p,1) + ...
          pnewAll(elem(:,2))*lambda(p,2) + ...
          pnewAll(elem(:,3))*lambda(p,3);
     pbar  = pbar + weight(p)*pt;
end
pbar = sum(pbar.*area)/sum(area);
pnewAll = pnewAll - pbar;

%% Output
info = struct('itStep', ite ,'resu',resu(1:ite), 'resp',resp(1:ite), 'flag', flag, 'Vcycle_num', Vcycle_num, 'time', time);
% plot(1:ite, log10(resu(1:ite)), 1:ite, log10(resp(1:ite)));
u = unewAll;
p = pnewAll;
end

