


J=1
E=[0.1,0.1,0.1,0.1]

Eigenzustande=transpose([[1,1,0,0];[1,0,1,0];[1,0,0,1];[0,1,1,0];[0,1,0,1];[0,0,1,1]])


H=Hamilton(J,0,Eigenzustande)



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






[V,D]=eig(H);
V
D
