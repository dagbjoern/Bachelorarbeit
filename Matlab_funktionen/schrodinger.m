function dtx=schrodinger(t,x)
  global H_0;
  global w;
  global E;
  global Gitterkonstante;
  V=(Gitterkonstante/2)*E*diag([-cos(w*t)+sin(w*t),-cos(w*t)-sin(w*t) ,cos(w*t)-sin(w*t) , cos(w*t)+sin(w*t)]);
  dtx= -1j*(H_0+V)*x;
