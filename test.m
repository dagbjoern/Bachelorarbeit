
addpath('Matlab_funktionen')





#Eigenzustande=transpose([[1,1,0,0];[1,0,1,0];[1,0,0,1];[0,1,1,0];[0,1,0,1];[0,0,1,1]])



Sprungterme=1


Potential=[0.0,0.1,0.2,0.5,1.0]
Energien=[0.1 , 0.02 , 0.1 ,0.2 , 0.4 ,0.8]
b = 'cool'
Potential(1)
Gitterkonstante=1
Phasenverschiebung=0
Anzahl=1         #Anzahl der Perioden
Frequenz=[[0,1,1,2,1,2,3,4,5]
         ,[1,3,2,3,1,1,1,1,1]];

for i=1:length(Potential)
#  ['sweet' num2str(Potential(i)) 'cool' ]

  Energie_1=-Potential(i)
  Energie_2=Potential(i)
  Energie_3=-Potential(i)
  Energie_4=Potential(i)
  H_0=Hamilton_0(Sprungterme,[Energie_1,Energie_2,Energie_3,Energie_4]) % 1 für Energien in diagonale


#function H_f=H_F(H_0,E,phi,a,Anzahl,frequenz)
#function H_f=H_F(H_0,E,phi,r1,r2,Anzahl)
  for l = 1:length(Frequenz)
    for k = 1:length(Energien)
      Matrix=H_F(H_0,Energien(k),Phasenverschiebung,Gitterkonstante,Anzahl,(Frequenz(1,l)/Frequenz(2,l))*Potential(i));
      [V,D]=eig(Matrix)
      e=eig(Matrix);
      % if k==1
      %   eigenwerte=eval(['eigenwerte_E' num2str(k)]);
      % end
      % if k!=1
      %   eigenwerte=[eigenwerte,eval(['eigenwerte_E' num2str(k)])];
      % end
      real_V=real(V)
      imag_V=imag(V)
      save(['build/Eigenwerte_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*100) '_w=' num2str(Frequenz(1,l)) '%' num2str(Frequenz(2,l)) 'a.txt'],'e')
      save(['build/Realpart_Eigenvektoren_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*100) '_w=' num2str(Frequenz(1,l)) '%' num2str(Frequenz(2,l)) 'a.txt'],'real_V')
      save(['build/Imagpart_Eigenvektoren_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*100) '_w=' num2str(Frequenz(1,l)) '%' num2str(Frequenz(2,l)) 'a.txt'],'imag_V')
      Parameter=[Gitterkonstante,Anzahl,Phasenverschiebung];
      save(['build/Parameter_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*100) '_w=' num2str(Frequenz(1,l)) '%' num2str(Frequenz(2,l)) 'a.txt'],'Parameter')

      % assignin ('base',['eigenwerte_E' num2str(k)],e);
    end
  end
end


Energien=transpose(Energien*100);
Potential=transpose(Potential*100);
Frequenz=transpose(Frequenz);
save('build/Durchlaufende_Energien.txt','Energien')
save('build/Durchlaufende_Potentiale.txt','Potential')
save('build/Durchlaufende_Frequenzen.txt','Frequenz')

%
% abs(transpose(V(:,1))*V(:,1))^2
% abs(transpose(V(:,1))*V(:,2))^2
% abs(transpose(V(:,1))*V(:,3))^2
% abs(transpose(V(:,1))*V(:,4))^2
% abs(transpose(V(:,1))*V(:,5))^2

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
