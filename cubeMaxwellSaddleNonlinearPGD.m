%% CUBEMAXWELLSADDLE solves Maxwell type equations in a cube using lowest order element.
% This is a special case of div u = g being nozero.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%clear; close all;

%% Defacult setting
tic
[node,elem] = cubemesh([-1,1,-1,1,-1,1],2/2);
bdFlag = setboundary3(node,elem,'Dirichlet');
h = 2/2;
for k = 1:4
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    h = h/2;
end
%%
omega = 4;
pde.J1 = @(p) [sin(omega*p(:,1)).*cos(omega*p(:,2)).*sin(omega*p(:,3)), ...
              cos(omega*p(:,1)).*sin(omega*p(:,2)).*sin(omega*p(:,3)), ...
              2*cos(omega*p(:,1)).*cos(omega*p(:,2)).*cos(omega*p(:,3))]*omega^2;
pde.exactu = @(p)[0*p(:,1), 0*p(:,2), ...
    cos(omega*p(:,1)).*cos(omega*p(:,2)).*cos(omega*p(:,3))];
pde.g_D = pde.exactu;
pde.curlu = @(p) [-cos(omega*p(:,1)).*sin(omega*p(:,2)).*cos(omega*p(:,3)), ...
              sin(omega*p(:,1)).*cos(omega*p(:,2)).*cos(omega*p(:,3)), 0*p(:,3)]*omega;
pde.g = @(p) -cos(omega*p(:,1)).*cos(omega*p(:,2)).*sin(omega*p(:,3))*omega;
%parameters used to define reluctivity
a0 = 10; %10
a1 = 70; %10*n  73.89055
a2 = 1;  %1
pde.nu = @(s) a0 + a1*exp(-a2*s);
l = @(p) ((cos(omega*p(:,1)).*sin(omega*p(:,2)).*cos(omega*p(:,3))).^2 +...
    (sin(omega*p(:,1)).*cos(omega*p(:,2)).*cos(omega*p(:,3))).^2)*omega^2;
f1 = @(p) (-cos(omega*p(:,1)).*sin(omega*p(:,2)).*cos(omega*p(:,3))).*...
    (sin(omega*p(:,1)).*sin(omega*p(:,2)).*cos(omega*p(:,3))) + ...
    (sin(omega*p(:,1)).*cos(omega*p(:,2)).*cos(omega*p(:,3))).*...
    (cos(omega*p(:,1)).*cos(omega*p(:,2)).*cos(omega*p(:,3)));
f2 = @(p) (-cos(omega*p(:,1)).*sin(omega*p(:,2)).*cos(omega*p(:,3))).*...
    (-cos(omega*p(:,1)).*cos(omega*p(:,2)).*cos(omega*p(:,3))) + ...
    (sin(omega*p(:,1)).*cos(omega*p(:,2)).*cos(omega*p(:,3))).*...
    (-sin(omega*p(:,1)).*sin(omega*p(:,2)).*cos(omega*p(:,3)));
f3 = @(p) (-cos(omega*p(:,1)).*sin(omega*p(:,2)).*cos(omega*p(:,3))).*...
    (cos(omega*p(:,1)).*sin(omega*p(:,2)).*sin(omega*p(:,3))) + ...
    (sin(omega*p(:,1)).*cos(omega*p(:,2)).*cos(omega*p(:,3))).*...
    (-sin(omega*p(:,1)).*cos(omega*p(:,2)).*sin(omega*p(:,3)));
pde.J = @(p) feval(pde.nu,l(p).^(1/2)).*feval(pde.J1,p) +...
    a1*(-a2)*exp(-a2*l(p).^(1/2)).*l(p).^(-1/2).*cross([f1(p),f2(p),f3(p)],pde.curlu(p))*omega^3;

%%

%% Parameters
maxIt = 1; 
N = zeros(maxIt,1); 
energyErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);
uIuhErr = zeros(maxIt,1);
stepsize = 1;
beta = 1;
NAA = 1; % every NAA steps to update the preconditioner IV and IQ;
Nstart = 20; % starting from the Nstart iteration, update every NAA steps
Nmg = 200; % the steps for HodgeLaplacian mg in every iteration
maxTPDIt = 1000;
tol = 10^(-3);

resultsPGD.para = [a0,a1,a2];
resultsPGD.stepsize = stepsize;
resultsPGD.NAA = NAA;
resultsPGD.Nstart = Nstart;
resultsPGD.Nmg = Nmg;
resultsPGD.meshsize = 2^(k+1);
resultsPGD.tol = tol;
resultsPGD.omega = omega;

%% preconditioner
IV = @(u) u;
IQ = @(p) p;

%% Finite Element Method        
for k = 1:maxIt   

    %[node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
 
    %%%%%%%%%%%%%%%%%%%%%%%
    N = size(node,1);
    NT = size(elem,1);
    curlE = pde.curlu;
    [lambda,w] = quadpts3(2);
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:) ...
            + lambda(p,4)*node(elem(:,4),:);
        if isnumeric(curlE) % a constant vector
            curlEp = repmat(curlE,NT,1);
        else % function handel
            curlEp = curlE(pxy);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%

    elemSort = sort(elem,2);
    [elem2dof,edge] = dof3edge(elemSort);
    %curluh = ones(size(elem,1),1);

    uI = edgeinterpolate(pde.exactu,node,edge);
    %uold = uI;
    uold = edgeinterpolate(@(p) ones(size(p)),node,edge);%ones(size(unew));
    [Dlambda,volume] = gradbasis3(node,elemSort);
    % curl phi = 2*Dlambda_i cross Dlambda_j;
    clear curlPhi
    curlPhi(:,:,6) = 2*cross(Dlambda(:,:,3),Dlambda(:,:,4),2);
    curlPhi(:,:,1) = 2*cross(Dlambda(:,:,1),Dlambda(:,:,2),2);
    curlPhi(:,:,2) = 2*cross(Dlambda(:,:,1),Dlambda(:,:,3),2);
    curlPhi(:,:,3) = 2*cross(Dlambda(:,:,1),Dlambda(:,:,4),2);
    curlPhi(:,:,4) = 2*cross(Dlambda(:,:,2),Dlambda(:,:,3),2);
    curlPhi(:,:,5) = 2*cross(Dlambda(:,:,2),Dlambda(:,:,4),2);
    curlEhp = zeros(NT,3);
    for kk = 1:6
        curlEhp = curlEhp + ...
            repmat(uold(elem2dof(:,kk)),1,3).*curlPhi(:,:,kk);
    end
    pde.mu = a0 + a1*exp(-a2*sum(curlEhp.^2,2).^(1/2));

    option = [];
    [eqn,info] = MaxwellsaddleMat(node,elem,bdFlag,pde,option);
    NE = size(eqn.edge,1);
    NT = size(elemSort,1);% Ndof = max(elem2dof(:)); %N = size(node,1);

    ii = eqn.StifMat.ii;
    jj = eqn.StifMat.jj;
    index = 0;
    sA = eqn.StifMat.sA;
    for i = 1:6
        for j = i:6
            sA(index+1:index+NT) = sA(index+1:index+NT).*pde.mu;
            index = index + NT;
        end
    end
    index = 0;
    sAp = eqn.StifMat.sA;
    for i = 1:6
        for j = i:6
            sAp(index+1:index+NT) = sAp(index+1:index+NT);
            index = index + NT;
        end
    end
    diagIdx = eqn.StifMat.diagIdx;
    upperIdx = eqn.StifMat.upperIdx;
    isFreeNode = eqn.freeNode;
    isFreeEdge = eqn.freeEdge;

    A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),NE,NE);
    AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),NE,NE);
    A0 = A + AU + AU';
    clear AU

    A = sparse(ii(diagIdx),jj(diagIdx),sAp(diagIdx),NE,NE);
    AU = sparse(ii(upperIdx),jj(upperIdx),sAp(upperIdx),NE,NE);
    Ap = A + AU + AU';
    

    ub = zeros(size(uI));
    ub(~isFreeEdge) = uI(~isFreeEdge);
    f0 = eqn.f0;
    f = f0 - A0*ub;
    f = f(isFreeEdge);
    A  = A0(isFreeEdge,isFreeEdge);
    Ap = Ap(isFreeEdge,isFreeEdge);
    K  = eqn.K(isFreeNode,isFreeNode);
    G0  = eqn.G;
    G  = eqn.G(isFreeEdge,isFreeNode);
    grad = eqn.grad(isFreeEdge,isFreeNode);
    M0 = eqn.M;
    M = M0(isFreeEdge,isFreeEdge);
    Minve = spdiags(diag(M0).^(-1),0,NE,NE);
    DMinve = Minve(isFreeEdge,isFreeEdge);
    %g0 = eqn.g(isFreeNode);
    pold = zeros(size(G,2),1);
    Gt = G';
    g = eqn.g - eqn.G'*ub;
    g0 = g(isFreeNode);

    Lp = Gt*DMinve*G;
    %Lp = Gt*G;
    L = ichol(Lp);
    Lp0 = G0'*Minve*G0;
    NodeDof = size(isFreeNode,1);
    bdidx = zeros(NodeDof,1);
    bdidx(~isFreeNode) = 1;
    Tbd = spdiags(bdidx,0,NodeDof,NodeDof);
    T = spdiags(1-bdidx,0,NodeDof,NodeDof);
    LpD = T*Lp0*T + Tbd;

    %% Multigrid solver for the initial guess
    option.isFreeEdge = isFreeEdge;
    option.isFreeNode = isFreeNode; 
    option.solvermaxit = Nmg;
    option.solver = 'CG'; % 'Vcycle', 'mg', diag'
    option.coarsegridsolver = 'direct';%'direct'; % 'direct'
    option.printlevel = 0;
    option.smoother = 'JAC';
    %option.smoother = 'GS';
    option.smoothingstep = 8; %8 for 10, 10 for 20, 12 for 30
    option.smoothingparameter = 0.5;
    option.smoothingratio = 3;
    option.tol = 10^(-9);
    bigAinitial = [M,G; Gt,sparse(size(G,2),size(G,2))];
    [uoldtmep,pnew,info] = mgMaxwellsaddle(Ap,G,f,g0,node,elem,bdFlag,M,grad,option,[]);
    uold = uI;
    %uold(isFreeEdge) = uoldtmep;    
    
    resi = A*uold(isFreeEdge) - f;
    residualErr1 = zeros(200,1);
    residualErr1(1) = norm(resi)/sqrt(length(resi));
    ite = 1;
    residualErr2 = residualErr1;
    residualErr3 = residualErr2;

    volume = abs(simplexvolume(node,elem)); % uniform refine in 3D is not orientation preserved
    volumemu = volume.*(pde.mu);
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
                        [volume;volume;volume;volume]/4,[max(elem(:)),1]);
    %DMinv = spdiags(1./Mvlump(isFreeNode),0,sum(isFreeNode),sum(isFreeNode));
    %idmat = eye(size(uold(isFreeEdge),1),size(uold(isFreeEdge),1));

    Mvlumpmu = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
                        [volumemu;volumemu;volumemu;volumemu]/4,[max(elem(:)),1]);
    Mmu = spdiags(Mvlumpmu(isFreeNode),0,sum(isFreeNode),sum(isFreeNode));
    DMinv = spdiags(Mvlumpmu(isFreeNode)./(Mvlump(isFreeNode).^2),0,sum(isFreeNode),sum(isFreeNode));

    NumVcycle = [];%info.itStep;
    NumVcycle_ite = [NumVcycle];
    
    while abs(residualErr1(ite))>tol*residualErr1(1) && ite < maxTPDIt
    %while max(abs(resi))>tol && ite < maxTPDIt

        volumemu = volume.*(pde.mu);
        Mvlumpmu = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
            [volumemu;volumemu;volumemu;volumemu]/4,[max(elem(:)),1]);
        DMinv = spdiags(Mvlumpmu(isFreeNode)./(Mvlump(isFreeNode).^2),0,sum(isFreeNode),sum(isFreeNode));

        if mod(ite,NAA) == 1 || ite <= Nstart
            AA = (A + G*DMinv*Gt);
        end

        unew_tmp0 = A*uold(isFreeEdge) - f + G*(DMinv*(Gt*uold(isFreeEdge) - g0));

        %% Direct solver
        % unew = uold(isFreeEdge) - stepsize*((A+beta*G*Gt)\(unew_tmp)); 

        %% Multigrid solver
        clear option
        option.isFreeEdge = isFreeEdge;
        option.solvermaxit = Nmg;
        option.solver = 'CG'; % 'Vcycle', 'mg', diag'
        option.coarsegridsolver = 'direct';%'direct'; % 'direct'
        option.printlevel = 0;
        option.smoother = 'GSJAC';
        %option.smoother = 'GS';
        option.smoothingstep = 8; %8 for 10, 10 for 20, 12 for 30
        option.smoothingparameter = 0.5;
        option.smoothingratio = 3;
        option.tol = 10^(-2); %10^(-2) not convergent for a1=130

        [unew_tmp1,info] = mgHodgeLapELarge(AA,unew_tmp0,node,elem,bdFlag,option);
        NumVcycle = NumVcycle +info.itStep;
        NumVcycle_ite = [NumVcycle_ite,info.itStep];
        
        %% mg, not sure
%         unew_tmp2 = Gt*unew_tmp1;
%         option.H1mgoption.solvermaxit = 10;
%         option.H1mgoption.printlevel = 0;
%         unew_tmp2 = zeros(NodeDof,1);
%         unew_tmp2(isFreeNode) = Gt*unew_tmp1;
%         [unew_tmp2,H1info] = mg(LpD,unew_tmp2,elem,option.H1mgoption,[]);
%         unew_tmp2 = Minve*G0*unew_tmp2;
%         unew_tmp2 = unew_tmp2(isFreeEdge);

        unew_tmp2 = Gt*unew_tmp1;
        [unew_tmp2,fl,re,it] = pcg(Lp,unew_tmp2,10^(-9),500,L,L');
        %unew_tmp2(isFreeNode) = 0;
        %unew_tmp2 = (Gt*G)\unew_tmp2;
        
        unew_tmp2 = DMinve*G*unew_tmp2;
        %unew_tmp2 = G*unew_tmp2;
        ProjGrad = unew_tmp1 - unew_tmp2;

        unew = uold(isFreeEdge) - stepsize*ProjGrad;

        unewAll = uI;   % boundary values
        unewAll(isFreeEdge) = unew;

        % compute curl phi
        curlEhp = zeros(NT,3);
        for kk = 1:6
            curlEhp = curlEhp + ...
                repmat(unewAll(elem2dof(:,kk)),1,3).*curlPhi(:,:,kk);
        end
        pde.mu = a0 + a1*exp(-a2*sum(curlEhp.^2,2).^(1/2));

        index = 0;
        sA = eqn.StifMat.sA;
        for i = 1:6
            for j = i:6
                sA(index+1:index+NT) = eqn.StifMat.sA(index+1:index+NT).*pde.mu;
                index = index + NT;
            end
        end
        A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),NE,NE);
        AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),NE,NE);
        A0 = A + AU + AU';
        A  = A0(isFreeEdge,isFreeEdge);
        clear AU

        f = f0 - A0*ub;
        f = f(isFreeEdge);

        uold = uI; 
        uold(isFreeEdge) = unew;

        ite = ite+1;
        resi = unew_tmp1;
        residualErr2(ite) = sqrt(sum(unew_tmp1.*unew_tmp0)/length(unew_tmp1));
        residualErr1(ite) = sqrt(sum(unew_tmp0.*unew_tmp0)/length(unew_tmp0));
        residualErr3(ite) = max(abs(unew_tmp0));
        disp(residualErr1(ite));
    end
    clear curlPhi;
    u = unewAll;
    disp(ite)
    
    timeTotal = toc
    resultsPGD.residualErr1 = residualErr1;
    resultsPGD.residualErr2 = residualErr2;
    resultsPGD.residualErr3 = residualErr3;
    resultsPGD.u = uold;
    resultsPGD.ite = ite;
    resultsPGD.timeTotal = timeTotal;
    resultsPGD.ite = ite;
    resultsPGD.NumVcycle = NumVcycle;
    resultsPGD.NumVcycle_ite = NumVcycle_ite;
    resultsPGD.stp = 'residualErr1';
    %save('resultsTPD','resultsTPD');

    fprintf('\n\n # of DoFs = %d \n',length(u));
    % compute error
    energyErr(k) = getHcurlerror3ND(node,elem,pde.curlu,u);
    L2Err(k) = getL2error3ND(node,elem,pde.exactu,u);
    uIuhErr(k) = sqrt((u-uI)'*A0*(u-uI));        
    fprintf('||curl(u-u_h)|| is %g \n',energyErr(k))
    N(k) = length(u);
    h(k) = 1./(size(node,1)^(1/3)-1);   
    resultsPGD.L2Err = L2Err(1);
    resultsPGD.energyErr = energyErr(1);
end