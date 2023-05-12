function Counterflow()
clc;
close all;
pipe_hot = GetPipe_hot();
pipe_cold = GetPipe_cold();

x = 0:0.01:pipe_hot.L;

[Tx_hot2, Tx_cold2]=CalcEuler3(173, 166.1, 50, 35, x, pipe_hot, pipe_cold);

plot (x,Tx_hot2, x,Tx_cold2)
legend ('Tx_hot2', 'Tx_cold2')
end


function lambdatr = lambda_altshul(Re, Eps)
lambdatr=0.11*(Eps+68/Re)^(1/4);
end


% метод Эйлера без трения
function [Tx_hot2, Tx_cold2]=CalcEuler3(Thot_in, Tcold_out, Ghot, Gcold, x, pipe_hot, pipe_cold)

S_hot= pi * pipe_hot.din^2 / 4;
v_hot=Ghot/(pipe_hot.density*S_hot);

S_cold = pi * pipe_cold.din^2 / 4;
v_cold=Gcold/(pipe_cold.density*S_cold);

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);

Tx_hot2(1)=Thot_in;
Tx_cold2(1)=Tcold_out;
n=length(x);
for i=1:n-1
    t_hot = Tx_hot2(i);
    t_cold = Tx_cold2(i);
    q_hot = -k * (t_hot-t_cold);
    q_cold= k * (t_hot-t_cold);

    P=pi*(pipe_hot.din+pipe_cold.din)/2;

    f_hot=(q_hot*P/(Ghot*pipe_hot.c));
    f_cold=(-q_cold*P/(Gcold*pipe_cold.c));
    
    dx = x(i+1) - x(i);
        Tx_hot2(i+1)=t_hot+dx*f_hot;
        Tx_cold2(i+1)=t_cold+dx*f_cold;
end 
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

pipe_hot.k=(2*pipe_hot.lambda)/(pipe_hot.douter*log(pipe_hot.douter/pipe_hot.din)); % коэффициент теплопередачи

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
pipe_cold.Eta=0.0000129;% динамическая вязкость
pipe_cold.Ny=pipe_cold.Eta/pipe_cold.density; %кинематическая вязкость

pipe_cold.k=(2*pipe_cold.lambda)/(pipe_cold.douter*log(pipe_cold.douter/pipe_cold.din)); % коэффициент теплопередачи

end
