function [xk,saida,exitflag]=filtro(STR,x0,rodaQuadprog)
%-----------------------------------------------------
% General SQP original filter algorithm
% minimize an objective function subject to equality and inequality
% constraints
% Step is computed by SQP algorithm
% Subrotines:
% internoreto.m -> compute the step by SQP algorithm
% testereto.m: test the trial point by the filter criterion
% restauracaoPQS.m: restoration procedure
% steihaug.m: minimize a quadratic function subject to a ball
% atualizareto.m -> filter update
%
% Inputs
% STR: structure composed by the array and matrix of the optimization problem
% x0: start point
% rodaQuadprog: flag 1-0 to use quadprog or the Torrealba Augmented Lagrangian technique
%-----------------------------------------------------

STR.tol_otim=1e-5;
STR.tol_viab=1e-5;
nIter_int = 0;
F  = [];

alfa = 0.1; %0.01
xk = x0;
n = length(xk);
k = 1;
options=optimset('Algorithm','interior-point-convex','TolX',1e-8,'Tolfun',1e-8,'TolCon',1e-8,'Display','off');
fig = 0;
Hsp = speye(n);
k_max =50;% STR.KmaxFiltro;

%=============
FF = [];
HHteste = [];
NNdc = [];
%===============
%% Inicializacao
% dados referentes as restricoes lineares
l = STR.lb;
u = STR.ub;
% para eliminar essas variaveis do calculo da complementaridade
ufin=find(u < 1e10); lfin = find(l > -1e10);  indfin=intersect(ufin,lfin);

x = xk;
%calculos % avalia funcao objetivo e restricoes
[f,g,h,hteste,cineq,ceq,Jineq,Jeq,Hf,STR]=calculos(x,STR);
fant = f;
H = Hf;
ng = 0;% numero de restricoes nao lineares para cada serie sintetica

delta = 1e6*min(norm(g),h);
meq = length(ceq);
mineq = length(cineq);
lambda_eq = sparse(meq,1);  % multiplicadores de lagrange iniciais das restricoes de igualdade
lambda_ineq = sparse(mineq,1); % multiplicadores de lagrange iniciais das restricoes de caixa
%passo = 0;
normL = 10;
complemen = 10;
ndc=10;
pare=0;

%% Loop
while  pare==0
    
    if k>=k_max
        exitflag = -10;
        break
    end
    
    corrente = [f-alfa*h (1-alfa)*h];
    Ftemp = F; % Temporary filter
    
    %===========================================================================================================
    % Calcula xk+1 (chama o algoritmo interno)
    %===========================================================================================================
    %internoreto
    [x,f,g,h,hteste,ceq,Jeq,cineq,Jineq,sucesso,~,lambda_eq,lambda_ineq,~,STR,nIter_int] = internoreto(xk,H,corrente,Ftemp,f,g,Jeq,Jineq,ceq,cineq,...
        lambda_eq,lambda_ineq,delta,STR,rodaQuadprog,nIter_int);
    %======================================================================
    % Atualiza a hessiana do modelo quadratico
    %======================================================================
    H = Hf;
    
    
    
    %%
    % Atualizacao do filtro
    if f>=fant
        F = atualizareto(F,corrente,fig,FF);
    end
    varf = abs(f-fant)/f;
    fant = f;
    passo = norm(x-xk,inf);
    xk = x;
    
    if hteste < STR.tol_viab
        %======================================================================
        % Medida de estacionaridade
        %======================================================================
        beq = Jeq*g; % GKP - ver figura pz_filtro2 - maio 2020
        if rodaQuadprog == 1
            [pz,~,ext1,output] = quadprog(Hsp,[],[],[],Jeq,beq,l-x+g,u-x+g,[],options);
            nIter_int = nIter_int + output.iterations;
        else
            [pz,ext1,ik,~] = LA_mochila(Hsp,[],Jeq',beq,l-x+g,u-x+g);
            nIter_int = nIter_int + ik;
        end
        
        if ext1==1
            dc=pz-g;
            ndc=norm(dc,inf);
        end
    end
    
    
    
    %======================================================================
    % KKT conditions
    %======================================================================
    Jac = [Jeq;Jineq];
    lambdak = [lambda_eq; lambda_ineq];
    normL = norm(g+Jac'*lambdak);
    complemen = norm(lambda_ineq(indfin).*(lambda_ineq(indfin)>1e-5).*cineq(indfin),inf);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Criterio de parada
    if (ndc < STR.tol_otim && hteste < STR.tol_viab) || (normL < STR.tol_otim && complemen < STR.tol_otim && hteste < STR.tol_viab)
        pare = 1;
    end
    %======================================================================
    fprintf('\n')
    fprintf('%4d %9d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n',k, round(f), h, hteste, ndc, normL, complemen, passo, delta);
    delta = max(1,sqrt(2)*passo);
    k = k+1;
    [nF,pF] = size(F);
    
    FF = [FF f];
    HHteste = [HHteste h];
    %    NNorml = [NNorml normL];
    NNdc = [NNdc ndc];
    %    CComplemen = [CComplemen complemen];
    
    if sucesso~=0
        xk1 = x;
        fxk1 = f;
        hxk1 = hteste;
        iteracoes = k;
        exitflag=0;
        saida =[k, f, h, hteste, ndc, normL, complemen, passo, delta, nIter_int];
        %         disp('======= INSUCESSO NO FILTRO ===============')
        %         fprintf('Exitflag = %2d  no filtro. \n',exitflag)
        %         fprintf('Erro do tipo = %2d  no filtro. \n',sucesso)
        %         disp('===========================================')
        return
    end
    xk1 = x;
    fxk1 = f;
    hxk1 = hteste;
    iteracoes = k;
    exitflag = 1;
    % cpu = toc;
    saida =[k, f, h, hteste, ndc, normL, complemen, passo, delta, nIter_int];
end




