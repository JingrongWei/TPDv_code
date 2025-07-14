%% CONVERGENCE OF FAST SOLVERS FOR MIXED FINITE ELEMENT METHOD (RT0-P0) FOR DARCY'S EQUATIONS
%
% This example is to show the convergence of fast solvers of mixed finite
% element (RT0-P0) approximation of the Darcy's equations.
%
% Reference 
%
% L. Chen. Multigrid Methods for Constrained Minimization Problems and
% Application to Saddle Point Problems. Submitted,2014.

close all
clear variables
%% Setting
h = 1/8;
[node,elem] = squaremesh([-1,1,-1,1],h); 
mesh = struct('node',node,'elem',elem);
option.mg.solvermaxIt = 80;
option.mg.solver = 'Vcycle';
option.mg.tol = 1e-1;
option.mg.printlevel = 0;
option.rateflag = 1;
option = mfemoption(option);

L0 = 3;
maxIt = 1;
maxN = 1e7;

% pde = Darcydata0;
pde = DFdata1;
bdFlag = setboundary(node,elem,'Neumann');
% mfemDarcy(mesh,pde,option);

%% Generate an initial mesh 
for k = 1:L0
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%     if isfield(pde,'K') && isnumeric(pde.K)
%         pde.K = repmat(pde.K,4,1); % prolongate piecwise K to the fine grid
%     end    
end

%% Initialize err and time array
% error
errpL2 = zeros(maxIt,1);   
errpIphL2 = zeros(maxIt,1); 
erruL2 = zeros(maxIt,1);   
% erruHdiv = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); 

% info
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
itStep = zeros(maxIt,1);
time = zeros(maxIt,1);
flag = zeros(maxIt,1);
Vcycle_num = zeros(maxIt,1);
%% Finite Element Method        
for k = 1:maxIt
    rng(4);
    tic
    % solve the equation
%         [p,u,eqn,info] =  DFP0P1TPD(node,elem,bdFlag,pde,option);
%     [p,u,eqn,info] =  DFP0P1TPD_NL(node,elem,bdFlag,pde,option);
%     [p,u,eqn,info] =  DFP0P1TPD_NL1(node,elem,bdFlag,pde,option);
%     [p,u,eqn,info] =  DFP0P1TPD_hessian(node,elem,bdFlag,pde,option);
     [p,u,eqn,info] =  DFP0P1TPD_Newton(node,elem,bdFlag,pde,option);
%      [p,u,eqn,info] =  DFP0P1TPD_FP(node,elem,bdFlag,pde,option);
%      [p,u,eqn,info] =  DFP0P1TPD_NL2(node,elem,bdFlag,pde,option);
   time(k) = toc;
    u = reshape(u, size(u,1)/2, 2);
    
    % compute error
    if isfield(pde,'exactu') && isfield(pde,'f')
        erruL2(k) = getL2error(node,elem,pde.exactu,u);
%         erruHdiv(k) = getHdiverrorRT0(node,elem,pde.exactudiv,u,[]);
        uI = Lagrangeinterpolate(pde.exactu,node,elem,'P0');
        uIuh_diff = sqrt(sum((uI-u).^2, 2));
        area = simplexvolume(node,elem);
        erruIuh(k)=sqrt(dot((uIuh_diff).^2,area));
    end
    if isfield(pde,'exactp')
        errpL2(k) = getL2error(node,elem,pde.exactp,p);
        % interpolation
        pI = Lagrangeinterpolate(pde.exactp,node,elem,'P1');        
        pIph_diff = (pI(elem(:, 1)) - p(elem(:, 1)) + pI(elem(:, 2)) - p(elem(:, 2))...
            + pI(elem(:, 3)) - p(elem(:, 3)))/3;
        errpIphL2(k) = sqrt(dot((pIph_diff).^2,area));
    end
    
    % record solver information
    itStep(k) = info.itStep;
    flag(k) = info.flag;
    Vcycle_num(k) = info. Vcycle_num;

    % plot 
    N(k) = length(p) + length(u);
    h(k) = 2./(sqrt(size(node,1))-1);
%     if option.plotflag && N(k) < 2e3 % show mesh and solution for small size
%        figure(1);  
%        showresult(node,elem,p);    
%     end
    if N(k) > maxN
        break;
    end
   
    % refine mesh
    [node,elem, bdFlag] = uniformrefine(node,elem,bdFlag);
%     if isfield(pde,'K') && isnumeric(pde.K)
%         pde.K = repmat(pde.K,4,1); % prolongate to the fine grid
%     end
end
%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrateh2(h(1:k),errpIphL2(1:k),1,'-*','||p_I-p_h||_{\infty}',...
               h(1:k),errpL2(1:k),1,'k-+','||p-p_h||');
    subplot(1,2,2)
    showrateh2(h(1:k),erruL2(1:k),1,'k-*','|| u - u_h||',...
               h(1:k),erruIuh(1:k),1,'m-+','|| u_I - u_h||');
%                h(1:k),erruHdiv(1:k),1,'-+','||div(u - u_h)||',... 
end

%% Output
err = struct('h',h(1:k),'N',N(1:k),'pL2',errpL2(1:k),'pIphL2',errpIphL2(1:k),...
             'uL2',erruL2(1:k), 'uIuL2',erruIuh(1:k));
            
%% Display error and CPU time
disp('Table: Error')
colname = {'#Dof','h','||p-p_h||','||p_I-p_h||','||u-u_h||'};
disptable(colname,err.N,[],err.h,'%0.2e',err.pL2,'%0.5e',err.pIphL2,'%0.5e',...
                     err.uL2,'%0.5e');

%% plot time growth
N = [ 197121  787457  3147777];
It_TPD_old = [ 11 54 252];
It_TPD = [4 18 85];
It_TPDscale = [3 14 70];
It_TPDNL = [3.7 16 73];
It_TPDIN = [2.2 10 40];
It_Newton = [2.6 12 54];

r_TPD = showrate(N,It_TPD,1,'-*');
hold on
r_TPDs = showrate(N,It_TPDscale,1,'-*');
r_TPDNL = showrate(N,It_TPDNL,1,'k-x');
r_TPDIN = showrate(N,It_TPDIN,1,'g-x');
r_Newton =  showrate(N,It_Newton,1,'g-x');
% title(['Rate of time growth'],'FontSize', 154);
h_legend = legend('TPDv',['N^{' num2str(r_TPD,2) '}'],...
                  'TPD-$\mu$',['N^{' num2str(r_TPDs,2) '}'],...
                  'TPDv-IMEX',['N^{' num2str(r_TPDNL,2) '}'],...
                  'TPDv-IN',['N^{' num2str(r_TPDIN,2) '}'],...
                 'Newton',['N^{' num2str(r_Newton,2) '}'],...
                  'LOCATION','Best');
set(h_legend,'FontSize',16);
xlabel('N','FontSize',16);
ylabel('CPU time','FontSize',16);



