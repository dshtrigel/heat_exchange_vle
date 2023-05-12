function dT = CalcEulerTempHead(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

Chot=pipe_hot.c*Ghot;
Ccold=pipe_cold.c*Gcold;

m = 1/Chot+1/Ccold;
P = pi*(pipe_hot.din+pipe_cold.din)/2;

dT(1)=Thot_in - Tcold_in;
n=length(x);

for i=1:n-1
   t = dT(i);
   k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);
   q = -k*m*t;
   f=q*P;
   dx = x(i+1) - x(i);
   dT(i+1)=t+dx*f;
end
end 