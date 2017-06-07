function Unter=untermatrix(Matrix,position_m,position_n,grosse)
#    [m,n]=size(Matrix);
  #  quotienten=m/grosse
     Unter=Matrix((position_m-1)*grosse+1 : position_m*grosse, (position_n-1)*grosse+1 : position_n*grosse);
end
