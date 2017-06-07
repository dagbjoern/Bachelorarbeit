function H_f=H_F(H_0,E,phi,a,Anzahl,frequenz)
  %[V,D]=eig(H_0);
  [zeilen,spalten]=size(H_0);
  H_f=zeros((Anzahl*2+1)*spalten,(Anzahl*2+1)*spalten);
  size(H_f)
  for n=(-Anzahl):Anzahl
    	n;
    for m=(-Anzahl):Anzahl
      H_mn=H_m_n(H_0,n,m,E,phi,a);
      for k=1:spalten
        k;
        for j=1:spalten
          #j,m,k, n
          Matrixelement=H_m_n(H_0,n,m,E,phi,a);
          H_f((m+Anzahl)*spalten+j,(n+Anzahl)*(spalten)+k)=Matrixelement(j,k);
          if(j==k && n==m)
            H_f((m+Anzahl)*spalten+j,(n+Anzahl)*(spalten)+k)+=n*frequenz;
          end
        end
      end
    end
  end
  H_f;
end
#function H_mn=H_m_n(H_0,n,m,E,phi,a)
