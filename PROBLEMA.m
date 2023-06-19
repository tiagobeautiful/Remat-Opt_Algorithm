% Objective function and its derivatives

function [f] = PROBLEMA(x,STR,flag)

a = STR.a;
P = STR.P;

if flag == 0
    f = (1/2)*x'*P*x - a'*x;
elseif flag == 1
    f = P*x - a;
elseif flag == 2
    f = P;
else
    error('fudeu')
end