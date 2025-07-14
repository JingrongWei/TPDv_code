function [p,u,eqn,info] = DFP0P1TPD_Newton(node,elem,bdFlag,pde,option)


if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

N = size(node,1); 

maxIt = 1e3;
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
% uoldAll = ones(Nu,1);
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
M = sparse(1:Nu, 1:Nu, [K.*area;K.*area], Nu, Nu);

eqn = struct('M',M,'B',B,...
             'f',F(1:Nu),'g',F(Nu+1:end),...
             'Bff',B(freep,:),'Bfft',B(freep,:)', 'gf',F(Nu+freep));


Jinv = pde.Juinv(uoldAll);
MJinv = [1./area; 1./area].*Jinv;
S = eqn.Bff*( MJinv*(eqn.Bfft));
bigS = [S sparse(Np-1,1); sparse(1,Np-1) 1];

% set up multilevel structure for mg solver
setupOption.solver = 'NO';
setupOption.freeDof =  [true(Np-1, 1); false];
bigS = [S sparse(Np-1,1); sparse(1,Np-1) 1];
[~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(bigS,poldAll,elem,setupOption);

% A = [M (B(1:Np-1, :))'; B(1:Np-1, :) sparse(Np-1, Np-1)]; % direct solver
% bigu = A\F(1:end-1);
% unewAll = bigu(1:Nu);
% pnewAll = bigu(Nu+1:end);
% pnewAll(Np) = 0;


%% Error
resu = zeros(maxIt,1);
resp = zeros(maxIt,1);

%% Solve the saddle point problem
stepsize = 1;
Vcycle_num = 0;
Bp =  eqn.Bfft*poldAll(freep);
% h =        1/1, 1/2, 1/4, 1/8
% stepsize = 0.8, 0.3, 0.3, 0.3
tic
for ite = 1:maxIt
%     fprintf('Iteration %i \n', ite);
     % Newton's method
    
     ru = M*uoldAll + Bp - eqn.f;
%      rp = eqn.Bff*uoldAll - eqn.gf;
     rp = -(eqn.Bff*uoldAll - eqn.gf);
%       J = pde.Ju(uoldAll);
%      MJ = [area; area].*J;

      Jinv = pde.Juinv(uoldAll);
     MJinv = [1./area; 1./area].*Jinv;

%      bigr = [ru; rp];
% 
%      bigJ = [MJ eqn.Bff'; -eqn.Bff sparse(Np-1, Np-1)]; % exact Jacobian

%      bigdu = minres(bigJ,-bigr, [], 200);
%      bigdu = -bigJ\bigr;
% 
%      du = bigdu(1:Nu);
%      dp = bigdu(Nu+1:end);
% 
%      unewAll = uoldAll + stepsize*du;
%      pnewAll =  poldAll;   % boundary values
%      pnewAll(freep) =  poldAll(freep) + stepsize*dp;

     % Schur complement
     S = eqn.Bff*( MJinv*(eqn.Bfft));
     rp =  eqn.Bff*( MJinv*ru) + rp;

%      % multigrid
     [Ai, Bi, BBi] = mgmat(S, Ai, Res, Pro);
      bigS = [S sparse(Np-1,1); sparse(1,Np-1) 1];
     [dpAll, info] = mg(bigS, [rp;0], elem, option.mg,Ai,Bi,BBi,Res,Pro,isFreeDof);
     Vcycle_num  = Vcycle_num  + info.itStep;

     du = MJinv*(ru-eqn.Bfft*dpAll(freep));
     unewAll = uoldAll - stepsize*du;
     pnewAll =  poldAll;   % boundary values
     pnewAll(freep) =  poldAll(freep) - stepsize*dpAll(freep);


     %      dp = S\rp;
%      du = MJ\(ru-eqn.Bff'*dp);
%       unewAll = uoldAll - stepsize*du;
%      pnewAll =  poldAll;   % boundary values
%      pnewAll(freep) =  poldAll(freep) - stepsize*dp;    


   % M. Mass matrix for P0 element
    if ~isempty(pde.K) && ~isnumeric(pde.K)       % d is a function   
        unorm = sqrt(unewAll(1:NT).^2 + unewAll(NT+1:end).^2);
        K = pde.K(unorm); % K is a scalar
    end
    M = sparse(1:Nu, 1:Nu, [K.*area;K.*area], Nu, Nu);
    eqn.M = M;

    % compute error
    Bp =  eqn.Bfft*pnewAll(freep) ;
    resu(ite) = norm(M*unewAll + Bp - eqn.f);
    resp(ite) = norm(eqn.Bff*unewAll  - eqn.gf);
     
    if ite> 1 && max(resu(ite)/normf, resp(ite)/normg) < tol
        flag = 1;
        break;
    end

    uoldAll = unewAll;
    poldAll = pnewAll;

%% 
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

