function ThermalCalc()
clc;
close all;
pipe = GetPipe();

x = 0:200:pipe.L;

% Основной расчет
Tx = CalcShukhov2(120, 291, 35, x, pipe);
Tx2 = CalcEuler(120, 291, 35, x, pipe);


plot(x, Tx, x, Tx2);
end

function lambdatr = lambda_altshul(Re, Eps)

lambdatr=0.11*(Eps+68/Re)^(1/4);

end

function Tx =CalcEuler(Tfluid_in, Touter, G, x, pipe)

S = pi * pipe.din^2 / 4;
v=G/(pipe.density*S);
Eps=0.015/pipe.din;

Re=v*pipe.din/pipe.Ny;
lambdatr=lambda_altshul(Re, Eps);

Tx(1)=Tfluid_in;
n=length(x);

for i=1:n-1
    t = Tx(i);
    q = pipe.k * (t-Touter);
    f=(-(4*q)/(pipe.din*pipe.density*pipe.c*v))+((lambdatr*v)/(2*pipe.din*pipe.c));
    dx = x(i+1) - x(i);

    Tx(i+1)=t+dx*f;
  
end
end

% Выполнить расчет по формуле Шухова (без трения)
function Tx = CalcShukhov(Tfluid_in, Touter, G, x, pipe)

a=(pipe.k*pi*pipe.douter*pipe.L)/(G*pipe.c);
Tx=Touter+(Tfluid_in-Touter)*exp(-a*x/pipe.L);

end

%По Лурье с учетом трения
function Tx = CalcShukhov2(Tfluid_in, Touter, G, x, pipe)

S = pi * pipe.din^2 / 4;
v=G/(pipe.density*S);
Eps=0.015/pipe.din;

Re=v*pipe.din/pipe.Ny;
lambdatr=lambda_altshul(Re, Eps);
u=(lambdatr*pipe.density*v^3)/(8*pipe.k);
a=(pi*pipe.k*pipe.din)/(G*pipe.c);

Tx=Touter+u+(Tfluid_in-Touter-u)*exp(-a*x);


end

% Сформировать параметры трубы
function pipe = GetPipe()
pipe.diz=0.616; % внешний диаметр изоляции
pipe.L=4000; % длина

pipe.c=2483; % теплоемкость
pipe.din=0.498; % внутренний диаметр
pipe.douter=0.508; % внешний диаметр

pipe.lambda=0.00499; % коэффициент теплопроводности

pipe.k=(2*pipe.lambda)/(pipe.douter*log(pipe.douter/pipe.din)); % коэффициент теплопередачи

pipe.density = 419; % плотность спг
pipe.Eta=0.00002;% динамическая вязкость
pipe.Ny=pipe.Eta/pipe.density; %кинематическая вязкость

end