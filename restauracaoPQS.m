function [x,proib]=restauracaoPQS(ceq,cineq,Jeq,Jineq,xk,corrente,STR)

Frest = [];
eprest = 1;
gama0 = 0.0625;
gama1 = 0.25;
gama2 = 2;
eta1 = 0.1;
eta2 = 0.9;
krest = 0;
proib = 0;
kmax = 40;


m = length([ceq;cineq]); % numero total de restricoes
meq = length(ceq); % numero de restricoes de igualdade
mineq=m-meq; % numero de restricoes de desigualdade
gamatheta=min(0.001,1/(2*sqrt(m)));

ctemp = cineq;
Jtemp = Jineq;
ii = (ctemp <= 0);
ctemp(ii) = [];
Jtemp(ii,:) = [];
thetabar = [ceq;ctemp];
Jbar = [Jeq;Jtemp];
thetatemp = max(0,cineq); % inviabilidade para desigualdade
theta = [ceq;thetatemp];
J = [Jeq;Jineq];
nablaf = J'*theta; % gradiente da funcao objetivo

deltarest = max(100,norm(theta));

if isempty(thetabar)
    return
end

x = xk;
while norm(theta)>=0.95*corrente(2) && krest<kmax
    
    
    d = steihaug(Jbar'*thetabar,Jbar'*Jbar,deltarest); % minimiza o modelo quadrático na bola de raio delta
    
    x = xk+d;
    
    %calculos  % calcula as restricoes e jacobianas em x
    [~,~,~,~,cineq,ceq,Jineq,Jeq]=calculos(x,STR);
    
    ctemp = cineq;
    Jtemp = Jineq;
    ii = (ctemp <= 0);
    ctemp(ii) = [];
    Jtemp(ii,:) = [];
    thetabarmais = [ceq;ctemp];
    Jbarmais = [Jeq;Jtemp];
    thetatemp = max(0,cineq);
    thetamais = [ceq;thetatemp];
    Jmais = [Jeq;Jineq];
    
    
    if isempty(thetabarmais)
        theta = thetamais;
        J = Jmais;
        nablaf = J'*theta;
        return
    end
    
    ared = 0.5*norm(thetabar)^2-0.5*norm(thetabarmais)^2;
    pred = -thetabar'*Jbar*d-0.5*d'*(Jbar'*Jbar)*d;
    rho = ared/pred;
    
    %===========================================================
    % Teste do filtro
    %===========================================================
    if isempty(Frest)
        aceito = 1;
        n = 0;
    else
        n = length(Frest(:,1));
        t = 0;
        for i = 1:n
            j = 1;
            while j<=m && abs(thetamais(j))-max(0,abs(Frest(i,j))-gamatheta*norm(Frest(i,:)))>=eprest
                j = j+1;
            end
            if j<=m
                t = t+1;
            end
        end
        if t == n
            aceito = 1;
        else
            aceito=0;
        end
    end
    %=============================================================
    if aceito==1
        xk = x;
        theta = thetamais;
        thetabar = thetabarmais;
        Jbar = Jbarmais;
        J = Jmais;
        nablaf = J'*theta;
        
        if rho<eta1
            %==============================================================
            %    Atualiza F
            %==============================================================
            i=1;
            while i<=n
                j = 0;
                A1 = Frest(i,:);
                for s = 1:m
                    if abs(A1(s))-abs(theta(s))>=eprest
                        j = j+1;
                    end
                end
                if j==m
                    Frest(i,:) = []; %elimina os pontos dominados
                    n = n-1;
                else
                    i = i+1;
                end
            end
            Frest = [Frest;theta'];
            %==============================================================
        end
    else
        if rho>=eta1 %&& norm(d)<=delta
            xk = x;
            theta = thetamais;
            thetabar = thetabarmais;
            Jbar = Jbarmais;
            J = Jmais;
            nablaf = J'*theta;
        end
    end
    
    %==============================================================
    %    Atualiza delta
    %==============================================================
    if rho<eta1
        deltarest = 0.5*(gama0+gama1)*deltarest;
    elseif rho>=eta1 && rho<eta2
        deltarest = 0.5*(1+gama1)*deltarest;
    elseif rho>=eta2
        deltarest = 0.5*(1+gama2)*deltarest;
    end
    
    krest = krest+1;
    
end

