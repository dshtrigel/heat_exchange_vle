function dT = CalcAnalytTempHead(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)
Chot=pipe_hot.c*Ghot;
Ccold=pipe_cold.c*Gcold;
m = 1/Chot+1/Ccold;
F=pi*x*(pipe_hot.din+pipe_cold.din)/2;
k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);
dT=(Thot_in - Tcold_in)*exp(-k*m*F);
end 