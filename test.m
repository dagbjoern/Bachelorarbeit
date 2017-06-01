


J=1


#Eigenzustande=transpose([[1,1,0,0];[1,0,1,0];[1,0,0,1];[0,1,1,0];[0,1,0,1];[0,0,1,1]])


H_0=Hamilton_0(J)



zeros(size(H_0))
%function H_mn=H_m_n(H_0,n,m,E,phi,r1,r2)

H_m_n(H_0,3,2,2,0,4,2)

% function  vektor=CTC(i,j,zustand)  #Sprungterme
%       vektor=zustand;
%       if (zustand(i)==1  && zustand(j)==0 )
%         vektor(i)=0;
%         vektor(j)=1;
%       end
% end
#[V,D]=eig(H)
#V*D*V^-1
#V(:,3)
#Normieren(V(:,3))


[V,D]=eig(H_0);
V
D

transpose(V(:,4))*H_0*V(:,4)
