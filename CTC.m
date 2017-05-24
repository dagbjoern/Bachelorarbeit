function  vektor=CTC(i,j,zustand)  #Sprungterme
      vektor=zustand;
      if (zustand(i)==1  && zustand(j)==0 )
        vektor(i)=0;
        vektor(j)=1;
      end
end
