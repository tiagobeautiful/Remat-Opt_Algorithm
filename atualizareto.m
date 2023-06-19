% This function updates the filter, including the pair (f, h).
% The pairs in the filter F0 are the vertices of the forbidden region.
% f = f0 - alpha * h0 and h = (1 - alpha) * h0

function [F]=atualizareto(F0,corrente,fig,FF)
ep=1e-6;
if isempty(F0)
    F=corrente;
else
f=corrente(1);
h=corrente(2);
oti=F0(:,1); 
inf=F0(:,2); 
nf=length(oti); 
poti=oti(min(4,nf):end);
pinf=inf(min(4,nf):end);
nf1=length(oti);                                      %achtung   
if fig==1
    figure(20)                                      %achtung
    clf
    %clf(figure(11))
    hold on
    grid on
    maxFF=max(FF);
    abmin = min([poti;f]); abmax = max([poti;f;maxFF]);
    ormin = 0; ormax = max([pinf;h]);
    axis1 = [abmin-0.1*(abmax+1-abmin) abmax + 0.1*(abmax+1-abmin)];  
    axis2 = [-0.1*ormax ormax + 0.1*ormax];
    abmin = axis1(1); abmax = axis1(2);
    ormin = axis2(1); ormax = axis2(2);
    axis([axis1 axis2])
    title('Filter')
    xlabel('objective')
    ylabel('infeasibility')
    plot(oti,inf,'or');

    for j=1:nf1
       plot([oti(j) abmax],[inf(j) inf(j)],'--r')
       plot([oti(j) oti(j)],[inf(j) ormax],'--r')
       %pause(.3)                                              %achtung
    end
    plot(f,h,'om');
    plot([f abmax],[h h],'--m')
    plot([f f],[h ormax],'--m')
    %pause(0.7)                                                %achtung

    dom=0;
    j=nf;
    while j>0 && f-F0(j,1)<=ep
        j=j-1;
    end
    while (j+dom+1)<=nf && h-F0(j+dom+1,2)<=ep
        dom=dom+1;
    end
    nf=nf+1-dom;
    if nf>j+1
        if dom==0
            for i=nf:-1:j+2
                F(i,1)=F0(i-1,1);
                F(i,2)=F0(i-1,2);
            end
               % F=Ordena([F0;corrente]);
        elseif dom>=1
            for i=j+2:nf
                F(i,1)=F0(i+dom-1,1);
                F(i,2)=F0(i+dom-1,2);
            end
        end
     end
    for i=1:j
            F(i,1)=F0(i,1);
            F(i,2)=F0(i,2);
    end

    F(j+1,1)=f;
    F(j+1,2)=h;
    plot(F(:,1),F(:,2),'ob');
    for j=1:nf
       plot([F(j,1) abmax],[F(j,2) F(j,2)],'--b')
       plot([F(j,1) F(j,1)],[F(j,2) ormax],'--b')
       %pause(.3)
    end
else
    dom=0;
    j=nf;
    while j>0 && f-F0(j,1)<=ep
        j=j-1;
    end
    while (j+dom+1)<=nf && h-F0(j+dom+1,2)<=ep
        dom=dom+1;
    end
    nf=nf+1-dom;
    if nf>j+1
        if dom==0
            for i=nf:-1:j+2
                F(i,1)=F0(i-1,1);
                F(i,2)=F0(i-1,2);
            end
%            F=Ordena([F0;corrente]);
        elseif dom>=1
            for i=j+2:nf
                F(i,1)=F0(i+dom-1,1);
                F(i,2)=F0(i+dom-1,2);
            end
        end
     end
    for i=1:j
            F(i,1)=F0(i,1);
            F(i,2)=F0(i,2);
    end

    F(j+1,1)=f;
    F(j+1,2)=h;
end
end
end
