
#addpath('Matlab_funktionen')
addpath('C:\Users\daghe\Desktop\Uni\Bachelorarbeit\Matlab_funktionen')


function En(x)
  global E;
  E = x;
endfunction

function Ham(x)
  global H_0;
  H_0 = x;
endfunction

function Fr(x)
  global w;
  w = x;
endfunction


#Eigenzustande=transpose([[1,1,0,0];[1,0,1,0];[1,0,0,1];[0,1,1,0];[0,1,0,1];[0,0,1,1]])



global E
global w
global H_0
global hbar=1

Sprungterme=1


Potential=[1,5,10]
Energien=[1,5,10]
Frequenz=[0.5,1,2,5,10]


Anzahl=[10,25]      #Anzahl der Perioden
b = 'cool'
global Gitterkonstante=1
Phasenverschiebung=0
%zeitentwicklung
t_isode=linspace(0,20,400)

#Frequenz=linspace(0,4,1000)
#[t,x]=rk4('test_rkt',[0,1],0.5)
#figure(1)
#plot(t,x)
#print figure1.pdf
#Frequenz=round_nur_besser(Frequenz,3)
test=0
for i=1:length(Potential)
   Energie_1=-Potential(i);
   Energie_2=Potential(i);
   Energie_3=-Potential(i);
   Energie_4=Potential(i);
  Ham(Hamilton_0(Sprungterme,[Energie_1,Energie_2,Energie_3,Energie_4])); % 1 für Energien in diagonale
  H_0_e=eig(H_0);
  [H_0_V,D]=eig(H_0);
  save(['Parameter/eigenvektoren_von_H_0_fur_a=' num2str(Potential(i)*100) '.txt'],'H_0_V')
  save(['Parameter/eigenwerte_von_H_0_fur_a=' num2str(Potential(i)*100) '.txt'],'H_0_e')
  %
%   #function H_f=H_F(H_0,E,phi,a,Anzahl,frequenz)
%   #function H_f=H_F(H_0,E,phi,r1,r2,Anzahl)
  for l = 1:length(Frequenz)
    	     Frequenz(l)
         Fr(Frequenz(l))
        for k = 1:length(Energien)
          En(Energien(k))
          E
          w
          [t,x]=rk4('schrodinger',[0,20],H_0_V(:,1));
          x_lsode=lsode('schrodinger_Isode',[H_0_V(1,1),H_0_V(2,1),H_0_V(3,1),H_0_V(4,1),0,0,0,0],t_isode);
          real_x_lsode=x_lsode(:,1:4);
          imag_x_lsode=x_lsode(:,5:8);
          real_x=real(x);
          imag_x=imag(x);
          save(['build/Realpart_Eigenvektoren_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*10000) '_w=' num2str(Frequenz(l)*1000) '.txt'],'real_x')
          save(['build/Imagpart_Eigenvektoren_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*10000) '_w=' num2str(Frequenz(l)*1000) '.txt'],'imag_x')
          save(['build/Realpart_Eigenvektoren_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*10000) '_w=' num2str(Frequenz(l)*1000) 'lsode.txt'],'real_x_lsode')
          save(['build/Imagpart_Eigenvektoren_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*10000) '_w=' num2str(Frequenz(l)*1000) 'lsode.txt'],'imag_x_lsode')
          for j = 1:length(Anzahl)
            test=test+1
            i %/ length(Potential)
            l %/ length(Frequenz)
            k%/ length(Energien),
            j% /length(Anzahl)
            'fortschritt'
            test/(length(Potential)*length(Frequenz)*length(Energien)*length(Anzahl))
            Matrix=H_F(H_0,Energien(k),Phasenverschiebung,Gitterkonstante,Anzahl(j),Frequenz(l));
            [V,D]=eig(Matrix);
            e=eig(Matrix);
            epsilon=e((size(e)/2-1):(size(e)/2+2));
            V_epsilon=V(:,(size(e)/2-1):(size(e)/2+2));
            real_phi=real(V_epsilon);
            imag_phi=imag(V_epsilon);
            save(['build/epsilon_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*10000) '_w=' num2str(Frequenz(l)*10000) '_N=' num2str(Anzahl(j) ) '.txt'],'epsilon')
            save(['build/Realpart_Eigenzustande_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*10000) '_w=' num2str(Frequenz(l)*1000) '_N=' num2str(Anzahl(j) ) '.txt'],'real_phi')
            save(['build/Imagpart_Eigenzustande_fur_a=' num2str(Potential(i)*100) '_E=' num2str(Energien(k)*10000) '_w=' num2str(Frequenz(l)*1000) '_N=' num2str(Anzahl(j) ) '.txt'],'imag_phi')
            Parameter=[Gitterkonstante,Anzahl(j),Phasenverschiebung];
            save(['Parameter/Parameter_fur_a=' num2str(Potential(i)*100) '_w=' num2str(Frequenz(l)*1000) '_E=' num2str(Energien(k)*10000) '_N=' num2str(Anzahl(j)) '.txt'],'Parameter')
          end
        end
   end
 end

 save(['build/Zeit_RK4.txt'], 't')
 save(['build/Zeit_lsode.txt'],'t_isode')
%

Energien=transpose(Energien*10000);
Potential=transpose(Potential*100);
Frequenz=transpose(Frequenz*1000);
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
