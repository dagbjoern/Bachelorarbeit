
addpath('Matlab_funktionen')





#Eigenzustande=transpose([[1,1,0,0];[1,0,1,0];[1,0,0,1];[0,1,1,0];[0,1,0,1];[0,0,1,1]])


Energie_1=
Energie_2=0
Energie_3=0
Energie_4=0
Sprungterme=1
H_0=Hamilton_0(Sprungterme,[Energie_1,Energie_2,Energie_3,Energie_4]) % 1 für Energien in diagonale


#function H_f=H_F(H_0,E,phi,a,Anzahl,frequenz)
#function H_f=H_F(H_0,E,phi,r1,r2,Anzahl)
Gitterkonstante=1
Frequenz=5
Anzahl=1         #Anzahl der Perioden
Phasenverschiebung=0
Energien=[]

% Variation der Energien von 0-4
maximale_E=5
for k = 1:1:maximale_E
  Energie=0+(k-1)/maximale_E;
  Energien=[Energien;Energie]
  Matrix=H_F(H_0,Energie,Phasenverschiebung,Gitterkonstante,Anzahl,Frequenz)
  e=eig(Matrix);
    assignin ('base',['eigenwerte_E' num2str(k)],e);
  if k==1
    eigenwerte=eval(['eigenwerte_E' num2str(k)]);
  end
  if k!=1
  eigenwerte=[eigenwerte,eval(['eigenwerte_E' num2str(k)])];
  end
end
Parameter=[Gitterkonstante,Frequenz,Anzahl,Phasenverschiebung]

save('build/Eigenwerte_Energien.txt','eigenwerte')
save('build/Parameter.txt','Parameter')
save('build/Energien.txt','Energien')


[V,D]=eig(Matrix)

abs(transpose(V(:,1))*V(:,1))^2
abs(transpose(V(:,1))*V(:,2))^2
abs(transpose(V(:,1))*V(:,3))^2
abs(transpose(V(:,1))*V(:,4))^2
abs(transpose(V(:,1))*V(:,5))^2

% Variation der Schritte
% for k = 1:1:5
%   e=eig(H_F(H_0,1,0,a,k,1))
%     assignin ('base',['eigenwerte_E' num2str(k)],e);
%   if k==1
%     eigenwerte=eval(['eigenwerte_E' num2str(k)])
%   end
%   if k!=1
%   eigenwerte=[eigenwerte,eval(['eigenwerte_E' num2str(k)])]
%   end
% end
% save('Eigenwerte_Große.txt','eigenwerte')


% Variation der frequenz
% for k = 1:1:5
%   e=eig(H_F(H_0,1,0,a,1,k));
%     assignin ('base',['eigenwerte_E' num2str(k)],e);
%   if k==1
%     eigenwerte=eval(['eigenwerte_E' num2str(k)]);
%   end
%   if k!=1
%   eigenwerte=[eigenwerte,eval(['eigenwerte_E' num2str(k)])];
%   end
% end
% save('Eigenwerte_Frequenz.txt','eigenwerte')
%
%
% Variation der Gitterkonstante
% for k = 1:1:5
%   e=eig(H_F(H_0,1,0,k,1,1));
%     assignin ('base',['eigenwerte_E' num2str(k)],e);
%   if k==1
%     eigenwerte=eval(['eigenwerte_E' num2str(k)]);
%   end
%   if k!=1
%   eigenwerte=[eigenwerte,eval(['eigenwerte_E' num2str(k)])];
%   end
% end
% save('Eigenwerte_Gitter.txt','eigenwerte')
%
%
%
%
%
%
%
%
% #function Unter=untermatrix(Martix,position_m,position_n,grosse)
% % U_11=untermatrix(H_f,1,1,4)
% % U_22=untermatrix(H_f,2,2,4)
% % U_33=untermatrix(H_f,3,3,4)
% % U_12=untermatrix(H_f,1,2,4)
% % U_21=untermatrix(H_f,2,1,4)
%
%
%
%
% #H_m_n(H_0,3,2,2,0,a)
% % function  vektor=CTC(i,j,zustand)  #Sprungterme
% %       vektor=zustand;
% %       if (zustand(i)==1  && zustand(j)==0 )
% %         vektor(i)=0;
% %         vektor(j)=1;
% %       end
% % end
% % #[V,D]=eig(H)
% % #V*D*V^-1
% % #V(:,3)
% % #Normieren(V(:,3))
%
%
% % [V,D]=eig(H_0);
% % V
% % D
% %
% % transpose(V(:,4))*H_0*V(:,4)
