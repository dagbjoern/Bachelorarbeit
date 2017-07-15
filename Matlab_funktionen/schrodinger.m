function dtx=schrodinger(t,x)
  global H_0;
  global w;
  global E;
  global Gitterkonstante;
  global hbar=1
  V=(Gitterkonstante/2)*E*diag([-cos(w*t)+sin(w*t),-cos(w*t)-sin(w*t) ,cos(w*t)-sin(w*t) , cos(w*t)+sin(w*t)]);
  dtx= -(1j/hbar)*(H_0+V)*x;
