function H_mn=H_m_n(H_0,n,m,E,phi,r1,r2)
H_mn=zeros(size(H_0))
  if(n==m)
    H_mn=+H_0
  end  % function
  if(n+1==m)
    H_mn=diag(ones(size(H_0),1)*(E/2)*exp(1j*phi)*(r1+1j*r2))
  end
  if(n-1==m)
    H_mn=+diag(ones(size(H_0),1)*(E/2)*exp(-1j*phi)*(r1-1j*r2))
  end
end
