function wert= round_nur_besser(X,nachkommerstelle)
  Y=X*10^nachkommerstelle;
  wert=round(Y);
  wert=wert/10^nachkommerstelle;
end
