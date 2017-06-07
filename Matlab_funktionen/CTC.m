function  vektor=CTC(i,j,zustand)  #Sprungterme und Eigenzustande
    if (i==j)
      vektor=zustand
    end
    if (zustand(i)==1  && zustand(j)==0 )
       vektor(i)=0;
        vektor(j)=1;
      end
end
