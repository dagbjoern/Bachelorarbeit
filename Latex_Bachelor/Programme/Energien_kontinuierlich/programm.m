
#addpath('Matlab_funktionen')
addpath('C:\Users\daghe\Desktop\Uni\Bachelorarbeit\Matlab_funktionen')

mkdir Parameter
mkdir build




#Eigenzustande=transpose([[1,1,0,0];[1,0,1,0];[1,0,0,1];[0,1,1,0];[0,1,0,1];[0,0,1,1]])



Sprungterme=1
global hbar=1


#Energien=[0.1 , 0.02 , 0.1 ,0.2 , 0.4 ,0.8]
Potential=[0.5,1.0,2.0]
Energien=linspace(0,6,150)
Energien=round(Energien*10000)/10000
Anzahl=[1,3,6,10]       #Anzahl der Perioden


b = 'cool'
Gitterkonstante=1
Phasenverschiebung=0

#Frequenz=linspace(0,4,1000)
Frequenz=[1.0,2.0]

#Frequenz=round_nur_besser(Frequenz,3)
fortschritt=0

for i=1:length(Potential)
#  ['sweet' num2str(Potential(i)) 'cool' ]
  Energie_1=-Potential(i);
  Energie_2=Potential(i);
  Energie_3=-Potential(i);
  Energie_4=Potential(i);
  H_0=Hamilton_0(Sprungterme,[Energie_1,Energie_2,Energie_3,Energie_4]); % 1 für Energien in diagonale
  H_0_e=eig(H_0)
  save(['Parameter/eigenwerte_von_H_0_fur_a=' num2str(Potential(i)*10000) 'a.txt'],'H_0_e')
#function H_f=H_F(H_0,E,phi,a,Anzahl,frequenz)
#function H_f=H_F(H_0,E,phi,r1,r2,Anzahl)
    for l = 1:length(Frequenz)
      for j = 1:length(Anzahl)
        for k = 1:length(Energien)
          fortschritt=fortschritt+1;
          i %/ length(Potential)
          l %/ length(Frequenz)
          j% /length(Anzahl)
          k%/ length(Energien),
          'fortschritt'
          fortschritt/(length(Potential)*length(Frequenz)*length(Energien)*length(Anzahl))
          Matrix=H_F(H_0,Energien(k),Phasenverschiebung,Gitterkonstante,Anzahl(j),Frequenz(l));
          #[V,D]=eig(Matrix);
          #real_V=real(V);
          #imag_V=imag(V);
          e=eig(Matrix);
          assignin ('base',['eigenwerte_E' num2str(l)],e);
          if k==1
            eigenwerte=eval(['eigenwerte_E' num2str(l)]);
          end
          if k!=1
            eigenwerte=[eigenwerte,eval(['eigenwerte_E' num2str(l)])];
          end
          #save(['build/Realpart_Eigenvektoren_fur_a=' num2str(Potential(i)*10000) '_E=' num2str(Energien(k)*10000) '_w=' num2str(Frequenz(l)*10000) 'a.txt'],'real_V')
          #save(['build/Imagpart_Eigenvektoren_fur_a=' num2str(Potential(i)*10000) '_E=' num2str(Energien(k)*10000) '_w=' num2str(Frequenz(l)*10000) 'a.txt'],'imag_V')
          % assignin ('base',['eigenwerte_E' num2str(k)],e);
        end
        save(['build/Eigenwerte_fur_a=' num2str(Potential(i)*10000) '_w=' num2str(Frequenz(l)*10000) '_N=' num2str(Anzahl(j)) '.txt'],'eigenwerte')
      end
  end
end


Energien=transpose(Energien*10000);
Potential=transpose(Potential*10000);
Frequenz=transpose(Frequenz*10000);
save('build/Durchlaufende_Energien.txt','Energien')
save('build/Durchlaufende_Potentiale.txt','Potential')
save('build/Durchlaufende_Frequenzen.txt','Frequenz')
save('build/Durchlaufende_N.txt','Anzahl')

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
