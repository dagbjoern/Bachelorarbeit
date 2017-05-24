function V_n = Normieren(Eigenvektor)
  sum=Eigenvektor(1)^2
  for i=2:length(Eigenvektor)
        sum=sum+Eigenvektor(i)^2
  end
  normierung=1/sqrt(sum)
  V_n=normierung*Eigenvektor
end
