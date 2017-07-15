function H_mn=H_m_n(H_0,n,m,E,phi,a)
global hbar
H_mn=zeros(size(H_0));
  if(n==m)
    H_mn=+H_0;
  end
  if(n+1==m)
    Hilfarray=ones(size(H_0),1)*(E*a/4)*exp(1j*phi);
    Hilfarray(1)=Hilfarray(1)*(-1-1j);
    Hilfarray(2)=Hilfarray(2)*(-1+1j);
    Hilfarray(3)=Hilfarray(3)*(1+1j);
    Hilfarray(4)=Hilfarray(4)*(1-1j);
    H_mn=-diag(Hilfarray);
  end
  if(n-1==m)
    Hilfarray=ones(size(H_0),1)*(E*a/4)*exp(-1j*phi);
    Hilfarray(1)=Hilfarray(1)*(-1+1j);
    Hilfarray(2)=Hilfarray(2)*(-1-1j);
    Hilfarray(3)=Hilfarray(3)*(1-1j);
    Hilfarray(4)=Hilfarray(4)*(1+1j);
    H_mn=-diag(Hilfarray);
  end
  H_mn;
end
