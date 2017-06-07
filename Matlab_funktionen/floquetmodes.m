function Phi=floquetmodes(V,n,w,t,c)
  l,k=size(V)
  Phi=zeros[l,1]
  for i=1:k
    for m=-n:n
    Phi=V[:,i]*exp(1j*m*t)
  end
end
