function [f,g,h,hteste,cineq,ceq,Jineq,Jeq,Hf,STR]=calculos(x,STR)
% Evaluates the objective function, constraints, 
% and their derivatives at point x and calculates 
% the value of the infeasibility measure at that point.
%%
A = STR.b';
b = STR.c;
l = STR.lb;
u = STR.ub;

[f] = PROBLEMA(x,STR,0);
% Gradiente da Funcao objetivo
[g] = PROBLEMA(x,STR,1);
% Hessiana da funcao objetivo
[Hf] = PROBLEMA(x,STR,2);

ceq = A*x-b;
cineq = [x-u;-x+l];
                            
% jacobiana das resticoes de caixa para teste do filtro
Hsp = speye(length(x));
Jineq = [Hsp;-Hsp];

% jacobiana das restricoes de igualdade
Jeq = A; % Jacobiana das retricoes de igualdade

% calculos para medida de inviabiliadde
thetatemp = min(1e20,max(0,cineq));
theta = [ceq;thetatemp]; 
h = norm(theta);
hteste = norm(theta,inf);










    
  




    
   
