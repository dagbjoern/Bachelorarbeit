
addpath('Matlab_funktionen')

#Gitterkonstante
a=10


J=1


#Eigenzustande=transpose([[1,1,0,0];[1,0,1,0];[1,0,0,1];[0,1,1,0];[0,1,0,1];[0,0,1,1]])


H_0=Hamilton_0(J)


#function H_f=H_F(H_0,E,phi,a,Anzahl,frequenz)
#function H_f=H_F(H_0,E,phi,r1,r2,Anzahl)


for k = 1:1:5
  E=0+(k-1)
  e=eig(H_F(H_0,k,0,a,1,2))
    assignin ('base',['eigenwerte_E' num2str(k)],e);
  if k==1
    eigenwerte=eval(['eigenwerte_E' num2str(k)])
  end
  if k!=1
  eigenwerte=[eigenwerte,eval(['eigenwerte_E' num2str(k)])]
  end
end
save('eigenwerte.txt','eigenwerte')



#function Unter=untermatrix(Martix,position_m,position_n,grosse)
% U_11=untermatrix(H_f,1,1,4)
% U_22=untermatrix(H_f,2,2,4)
% U_33=untermatrix(H_f,3,3,4)
% U_12=untermatrix(H_f,1,2,4)
% U_21=untermatrix(H_f,2,1,4)




#H_m_n(H_0,3,2,2,0,a)

'cool'
% function  vektor=CTC(i,j,zustand)  #Sprungterme
%       vektor=zustand;
%       if (zustand(i)==1  && zustand(j)==0 )
%         vektor(i)=0;
%         vektor(j)=1;
%       end
% end
% #[V,D]=eig(H)
% #V*D*V^-1
% #V(:,3)
% #Normieren(V(:,3))


% [V,D]=eig(H_0);
% V
% D
%
% transpose(V(:,4))*H_0*V(:,4)
