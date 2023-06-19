%Internal algorithm for calculating the step using PQS (straight filter).

function [x,f,g,h,hteste,ceq,Jeq,cineq,Jineq,sucesso,inviavel,lambda_eq,lambda_ineq,delta,STR,nIter_int] = ...
    internoreto(xk,H,corrente,F,fk,gk,Jeq,Jineq,ceq,...
    cineq,lambda_eq,lambda_ineq,delta,STR,rodaQuadprog,nIter_int)

% Dados
inviavel=0;
options=optimset('Algorithm','interior-point-convex','TolX',1e-8,'Tolfun',1e-8,'TolCon',1e-8,'Display','off');

cp = 0.0001;
xi = 0.8;
eta = 0.01;
gama = 0.5;
n = length(xk);
Hsp = speye(n);
sucesso = 0;
Jeqk = Jeq;
Jineqk = Jineq;
cineqk = cineq;
meq = length(ceq);

%% Passo de viabilidade
%==========================================================================
if rodaQuadprog
    [nk,~,EXITFLAG, output] = quadprog(Hsp,[],[],[],Jeqk,-ceq,STR.lb-xk,STR.ub-xk,[],options);
    aux_nIter = output.iterations;
    clear output
else
    [nk,EXITFLAG,aux_nIter,~] = LA_mochila(Hsp,[],Jeqk',-ceq,STR.lb-xk,STR.ub-xk);
end

if EXITFLAG==1 || EXITFLAG==0
    x1 = xk+nk;
    ceq1 = [STR.b'*x1-STR.c];
    cineq1 = [x1-STR.ub;-x1+STR.lb]; %tranforma caixa em restricoes de desigualdade
    theta1 = [ceq1; max(0,cineq1)]; % function c+
    h1 = norm(theta1); % infeasibility measure
end
%==========================================================================
if EXITFLAG==-2
    inviavel=1;
    [x,proib]=restauracaoPQS(ceq,cineq,Jeq,Jineq,xk,corrente,STR);
    [f,g,h,hteste,cineq,ceq,Jineq,Jeq,~,STR]=calculos(x,STR);
    tentativo = [f h];
    proib = testereto(F,tentativo,corrente,0);
    if proib==1
        sucesso = 4;
        return
    end
else
    %%  Passo de otimalidade
    pare = 0;
    cont = 0;
    x = xk;
    while pare==0
        if norm(nk,inf)>xi*delta
            if EXITFLAG==1 && h1<corrente(2)
                x = x1;
                [f,g,h,hteste,cineq,ceq,Jineq,Jeq,~,STR]=calculos(x,STR);
                pare = 1;
            else
                xk = x;
                [x,proib]=restauracaoPQS(ceq,cineq,Jeq,Jineq,xk,corrente,STR);
                [f,g,h,hteste,cineq,ceq,Jineq,Jeq,~,STR]=calculos(x,STR);
                tentativo = [f h];
                proib = testereto(F,tentativo,corrente,0);
                if proib==1
                    sucesso = 4;
                    return
                end
                pare = 1;
            end
        else
            %==========================================================================
            % Calculo do passo de otimalidade
            %==========================================================================
            [t,~,EXITFLAG,output,lambda] = quadprog(H,gk+H*nk,[],[],Jeqk,sparse(meq,1),max(STR.lb-xk-nk,-nk-delta*ones(n,1)),...
                min(STR.ub-xk-nk,-nk+delta*ones(n,1)),[],options);
            aux_nIter = aux_nIter + output.iterations;
            clear output
            %==========================================================================
            
            d = nk+t;
            x = xk+d;
            
            [f,g,h,hteste,cineq,ceq,Jineq,Jeq,~,STR]=calculos(x,STR);  % avalia as funcoes para fazer teste do filro
            
            tentativo = [f h];
            proib = testereto(F,tentativo,corrente,0);
            
            ared = fk-f;
            pred = -gk'*d-0.5*d'*H*d;
            
            if proib==1 || (pred>=cp*(h^2) && ared<eta*pred)
                delta = gama*delta;
                pare = 0;
            else
                pare = 1;
                lambda_eq = lambda.eqlin;
                lambda_ineq = [lambda.upper;lambda.lower];
            end
        end
        cont = cont + 1;
    end
end


nIter_int = nIter_int + aux_nIter + cont;
