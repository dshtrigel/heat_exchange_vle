function temperature_head()
clc;
close all;
pipe_hot = GetPipe_hot();
pipe_cold = GetPipe_cold();


x = 0:0.1:pipe_hot.L;

dT = CalcEule(173, 153, 50, 35, x, pipe_hot, pipe_cold);
dT2 = CalcTemperature(173, 153, 50, 35, x, pipe_hot, pipe_cold);


plot (x,dT, x, dT2);




end

function dT = CalcEule(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

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

function dT2 = CalcTemperature(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

Chot=pipe_hot.c*Ghot;
Ccold=pipe_cold.c*Gcold;
m = 1/Chot+1/Ccold;
F=pi*x*(pipe_hot.din+pipe_cold.din)/2;
k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);
dT2=(Thot_in - Tcold_in)*exp(-k*m*F);


end 


% Сформировать параметры трубы 1
function pipe_hot = GetPipe_hot()
pipe_hot.L=50; % длина

pipe_hot.din=0.10; % внутренний диаметр
pipe_hot.douter=0.108; % внешний диаметр
pipe_hot.wall=0.004; % толщина стенки

pipe_hot.lambda=209; % коэффициент теплопроводности

pipe_hot.alpha=52250;% коэффициент теплоотдачи
pipe_hot.c=5477; % теплоемкость
pipe_hot.density = 419; % плотность СПГ
pipe_hot.Eta=0.0000049;% динамическая вязкость
pipe_hot.Ny=pipe_hot.Eta/pipe_hot.density; %кинематическая вязкость

%pipe_hot.k = 1 /(1/pipe_hot.alpha + 1/(2*pipe_hot.lambda*log10(pipe_hot.douter/pipe_hot.din)));

%pipe_hot.k=(2*pipe_hot.lambda)/(pipe_hot.douter*log10(pipe_hot.douter/pipe_hot.din)); % коэффициент теплопередачи

end

% Сформировать параметры трубы 2
function pipe_cold = GetPipe_cold()
pipe_cold.L=50; % длина

pipe_cold.din=0.08; % внутренний диаметр
pipe_cold.douter=0.088; % внешний диаметр
pipe_cold.wall=0.004; % толщина стенки

pipe_cold.lambda=209; % коэффициент теплопроводности

pipe_cold.alpha=52250; % коэффициент теплоотдачи
pipe_cold.c=3950; % теплоемкость
pipe_cold.density =1349; % плотность азота
pipe_cold.Eta=0.000286;% динамическая вязкость
pipe_cold.Ny=pipe_cold.Eta/pipe_cold.density; %кинематическая вязкость

%pipe_cold.k=(1)/(1/pipe_hot.alpha + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha ); % коэффициент теплопередачи

end




