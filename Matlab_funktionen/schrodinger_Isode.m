function dtx=schrodinger_Isode(x,t)
  global H_0;
  global w;
  global E;
  global Gitterkonstante;
  global hbar;
    V=(Gitterkonstante/2)*(-E)*diag([-cos(w*t)+sin(w*t),-cos(w*t)-sin(w*t) ,cos(w*t)-sin(w*t) , cos(w*t)+sin(w*t)]);
    y=x(1:4)+i*x(5:8);
    H=-(i/hbar)*(H_0+V);
    dtx(1:4)=real(H*y);
    dtx(5:8)=imag(H*y);
