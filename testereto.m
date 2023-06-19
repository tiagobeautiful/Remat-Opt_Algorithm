% This function tests whether a tentative point is forbidden or not by the current filter (original filter).
% Input:
% F = temporary filter
% The filter pairs are of the form (f - alpha * h, (1 - alpha) * h)
% tentative = tentative point (f(x+), h(x+))
% current = current point (f(xk) - alpha * h(xk), (1 - alpha) * h(xk))
% alpha = factor defining the margin
% fig = shows the graph if fig = 1
% Output:
% proib = 0 if the point is accepted by the filter
% proib = 1 if the point is forbidden

function proib=testereto(F,tentativo,corrente,fig)

ep=0;
if isempty(F)
    nf=0;
else
    nf=length(F(:,1));
end
%F=Ordena(F);
f0=corrente(1);
h0=corrente(2);
f=tentativo(1);
h=tentativo(2);
%=========================================================================
%Gráfico
%=========================================================================
if nf~=0 && fig==1
    oti=F(:,1);
    inf=F(:,2);
    clf(figure(1))
    hold on
    grid on
    abmin = min([oti;f]); abmax = max([oti;f]);
    ormin = 0; ormax = max([inf;h;3]);
    axis1 = [abmin-0.1*(abmax+1-abmin) abmax + 0.1*(abmax+1-abmin)] ;
    axis2 = [-0.1*ormax ormax + 0.1*ormax];
    abmin = axis1(1); abmax = axis1(2);
    ormin = axis2(1); ormax = axis2(2);
    axis([axis1 axis2])
    title('Filtro')
    xlabel('objetivo')
    ylabel('inviabilidade')
    plot(oti,inf,'or');
    
    for j=1:nf
        plot([oti(j) abmax],[inf(j) inf(j)],'--r')
        plot([oti(j) oti(j)],[inf(j) ormax],'--r')
    end
    pause
    plot(f0,h0,'om');
    plot(f,h,'ok');
    plot([f0 abmax],[h0 h0],'--m')
    plot([f0 f0],[h0 ormax],'--m')
    pause
    
    %=========================================================================
    %Teste
    %=========================================================================
    
    if f-f0>=ep && h>=h0
        proib=1;
    else
        j=nf;
        while j>0 && f-F(j,1)<ep
            j=j-1;
        end
        if j>0 && h-F(j,2)>=ep
            proib=1;
        else
            proib=0;
        end
    end
else
    if f-f0>=ep && h-h0>=ep
        proib=1;
    else
        j=nf;
        while j>0 && f-F(j,1)<ep
            j=j-1;
        end
        if j>0 && h-F(j,2)>=ep
            proib=1;
        else
            proib=0;
        end
    end
end
end
