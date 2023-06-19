%% Augmented Lagrangian for knapsack
%% min   0.5 x'Px - a'x
% s.t.  b'x = c
%       l <= x <= u
% where x\in R^n, P=identidade(n,n), b=ones(1,n), a=[2, 3,...,n]',
% l=zeros(1,n), u = Constante.
%%
function [x,sucesso,ik,t] = LA_mochila(P,a,b,c,l,u)

r=1;
beta = 1;
lamb=1;
kmax=6000;
sucesso = -2;
prec = 10^(-8);
n = length(u);
%lamb = (1/c)*sum(sqrt(g))^2;
xk0 = (l+u)./2;
erro2 = 1;
erro1 = 1;
tic;    % Time
somak = 0;

p = diag(P);
% b = p;
% y = P\a;%a./p;
z = b./p;
alfa=r/(1+r*(b'*z));
p1 = (beta*alfa*(b'*z) - 1)*z;
% p2 = y + (r*c - alfa*beta*b'*(y + r*c*z))*z;
p2 = (r*c - alfa*beta*b'*r*c*z)*z;

for k = 1 : kmax
    
    xb = p2 + lamb*p1;
    somak = somak + 1;
    
    %% Projecting xb on the box
    xk = max(l,min(xb,u));
    %% Atualizamos o multiplicador da restricao de igualdade
    rest= b'*xk-c;
    lamb = lamb + r*rest;
    
%     if abs(b'*xk)>0.9*abs(b'*xk0)
%         r=2*r;
%     end
        
    %% Criterio de parada
    erro1=norm(xk-xk0);
    erro2 = abs(rest);
    %pause
    if erro2 <= prec
        if erro1 <= prec
            sucesso = 1;
            break;
        end
    end
    xk0=xk;
end
t = toc();

x = xk;
ik = somak;
