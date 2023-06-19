function [d,j]=steihaug(g,B,delta)

% Minimization of a quadratic with ball constraint - Steihaug CG
% Nocedal/Wright p. 75

% Problem
%    min g'd + 0.5d'Bd
%    s.t |d|<=delta

% Data 
%    g: gradient 
%    B: symmetric matrix
%    delta: radius of the ball constraint

% Output
%    vector d, approximate solution of the problem

% Subroutine for 
%    qrlb: minimization of a quadratic with linear and ball constraints

% November 17, 2004 - Ademir and Elizabeth

% parameters
eps = 1e-4;   % precision
graph = 1;    % 1 for graphics, 0 otherwise
display = 1;  % 1 to display tables, 0 otherwise

% initialize
n = length(g); d = zeros(n,1); r = g; p =-r;
tol = eps*norm(r); jmax = n+1;
goon = 1; j = 1; 
if norm(r) <= eps
   goon=0;
end;

while j<jmax && goon
   curvature = p'*B*p;
   if curvature <=0
      % find t s.t. (d + tp) minimize the quadratic and |d + tp| = delta
      aa = p'*p;  bb = 2*d'*p; cc = (d'*d)-delta ^2;
      tt = roots([aa bb cc]);
      t = max(tt);
      d = d +t*p;
      goon = 0;
   else
     alpha = r'*r/curvature; dtemp = d; d = d + alpha*p;
     if norm(d)>=delta
        d = dtemp;
        aa = p'*p;  bb = 2*d'*p; cc = (d'*d)-delta^2;
        tt = roots([aa bb cc]);
        t=max(tt);
        d = d +t*p;
	goon = 0;
     else
        rtemp = r; r = r + alpha*B*p;
	if norm(r)<tol
	   goon = 0;
	end;
	beta = r'*r/(rtemp'*rtemp);
	p =-r + beta*p;
     end
   end
   j=j+1;
end


