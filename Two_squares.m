function Two_squares ()
clc;
close all;
pipe_hot = GetPipe_hot();
pipe_cold = GetPipe_cold();


x = 0:0.01:pipe_hot.L;

%[Tx_hot, Tx_cold]=CalcEuler(173, 153, 50, 35, x, pipe_hot, pipe_cold);
[Thot_out, Tcold_out]= FinalTemperatures(173, 153, 50, 35, x, pipe_hot, pipe_cold);
[Tx_hot2, Tx_cold2]=CalcEuler3(173, 153, 50, 35, x, pipe_hot, pipe_cold);
% [Thot_out2, Tcold_out2]= FinalTemperatures2(173, 153, 50, 35, x, pipe_hot, pipe_cold); 

plot (x,Thot_out, x,Tcold_out)
legend('Thot_out', 'Tcold_out')

figure
%plot(x,Tx_hot, 'LineWidth', 1)
%hold on
%plot(x,Tx_cold, 'LineWidth', 1)
%hold on
plot(x,Thot_out, 'LineWidth', 3)
hold on
plot(x,Tcold_out, 'LineWidth', 3)
hold on
plot (x,Tx_hot2, 'LineWidth', 1); 
hold on
plot (x,Tx_cold2, 'LineWidth', 1);
hold on
legend('Thot_out', 'Tcold_out', 'Tx_hot2', 'Tx_hot2')
hold off

end


function lambdatr = lambda_altshul(Re, Eps)

lambdatr=0.11*(Eps+68/Re)^(1/4);

end
% Метод Эйлера
function [Tx_hot, Tx_cold]=CalcEuler(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

S_hot= pi * pipe_hot.din^2 / 4;
v_hot=Ghot/(pipe_hot.density*S_hot);
Eps_hot=0.015/pipe_hot.din;
Re_hot=v_hot*pipe_hot.din/pipe_hot.Ny;
lambdatr_hot=lambda_altshul(Re_hot, Eps_hot);

S_cold = pi * pipe_cold.din^2 / 4;
v_cold=Gcold/(pipe_cold.density*S_cold);
Eps_cold=0.015/pipe_cold.din;
Re_cold=v_cold*pipe_cold.din/pipe_cold.Ny;
lambdatr_cold=lambda_altshul(Re_cold, Eps_cold);

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + 1/pipe_cold.alpha);

Tx_hot(1)=Thot_in;
Tx_cold(1)=Tcold_in;
n=length(x);
for i=1:n-1
    t_hot = Tx_hot(i);
    t_cold = Tx_cold(i);
   
    q_hot = -k * (t_hot-t_cold);
    q_cold= k * (t_hot-t_cold);

    P=pi*pipe_hot.douter;

    f_hot=((q_hot)*P/(Ghot*pipe_hot.c))+((lambdatr_hot*v_hot^2)/(2*pipe_hot.din*pipe_hot.c));
    f_cold=((q_cold)*P/(Gcold*pipe_cold.c))+((lambdatr_cold*v_cold^2)/(2*pipe_cold.din*pipe_cold.c));

%     f_hot=(-(4*q_hot)/(pipe_hot.din*pipe_hot.density*pipe_hot.c*v_hot))+((lambdatr_hot*v_hot^2)/(2*pipe_hot.din*pipe_hot.c));
%     f_cold=(-(4*q_cold)/(pipe_cold.din*pipe_cold.density*pipe_cold.c*v_cold))+((lambdatr_cold*v_cold^2)/(2*pipe_cold.din*pipe_cold.c));  
    
    dx = x(i+1) - x(i);
    Tx_hot(i+1)=t_hot+dx*f_hot;
    Tx_cold(i+1)=t_cold+dx*f_cold;
end 
end

% метод Эйлера без трения
function [Tx_hot2, Tx_cold2]=CalcEuler3(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

% S_hot= pi * pipe_hot.din^2 / 4;
% v_hot=Ghot/(pipe_hot.density*S_hot);
% 
% S_cold = pi * pipe_cold.din^2 / 4;
% v_cold=Gcold/(pipe_cold.density*S_cold);

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);

Tx_hot2(1)=Thot_in;
Tx_cold2(1)=Tcold_in;
n=length(x);
for i=1:n-1
    t_hot = Tx_hot2(i);
    t_cold = Tx_cold2(i);
    q_hot = -k * (t_hot-t_cold);
    q_cold= k * (t_hot-t_cold);

    P=2*pipe_hot.a+2*pipe_hot.wall;

    f_hot=((q_hot)*P/(Ghot*pipe_hot.c));
    f_cold=((q_cold)*P/(Gcold*pipe_cold.c));

%     f_hot=(-(4*q_hot)/(pipe_hot.din*pipe_hot.density*pipe_hot.c*v_hot));
%     f_cold=(-(4*q_cold)/(pipe_cold.din*pipe_cold.density*pipe_cold.c*v_cold));
    
    dx = x(i+1) - x(i);
    Tx_hot2(i+1)=t_hot+dx*f_hot;
    Tx_cold2(i+1)=t_cold+dx*f_cold;
end 
end


% вывод по Исаченко
function [Thot_out, Tcold_out]= FinalTemperatures(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);

Chot=pipe_hot.c*Ghot;
Ccold=pipe_cold.c*Gcold;
m = 1/Chot+1/Ccold;

% P=2*pipe_hot.a+2*pipe_hot.wall;
% F=x*P;

F=2*(x*pipe_hot.a+pipe_hot.a*2*pipe_hot.wall+x*pipe_hot.wall);

Thot_out=Thot_in-(Thot_in-Tcold_in)*(1-exp(-k*F*m))/(1+Chot/Ccold);
Tcold_out=Tcold_in+(Thot_in-Tcold_in)*(1-exp(-k*F*m))/(1+Ccold/Chot);

end

  


% Сформировать параметры трубы 1
function pipe_hot = GetPipe_hot()
pipe_hot.L=50; % длина

pipe_hot.a=0.10; % внутренний стенка

pipe_hot.wall=0.004; % толщина стенки

pipe_hot.lambda=209; % коэффициент теплопроводности

pipe_hot.alpha=52250;% коэффициент теплоотдачи
pipe_hot.c=5477; % теплоемкость
pipe_hot.density = 419; % плотность СПГ
pipe_hot.Eta=0.0000049;% динамическая вязкость
pipe_hot.Ny=pipe_hot.Eta/pipe_hot.density; %кинематическая вязкость

%pipe_hot.k=(2*pipe_hot.lambda)/(pipe_hot.douter*log(pipe_hot.douter/pipe_hot.din)); % коэффициент теплопередачи

end

% Сформировать параметры трубы 2
function pipe_cold = GetPipe_cold()
pipe_cold.L=50; % длина

pipe_cold.a=0.1; % внутренний диаметр

pipe_cold.wall=0.004; % толщина стенки

pipe_cold.lambda=209; % коэффициент теплопроводности

pipe_cold.alpha=52250; % коэффициент теплоотдачи
pipe_cold.c=3950; % теплоемкость
pipe_cold.density =1349; % плотность азота
pipe_cold.Eta=0.0000129;% динамическая вязкость
pipe_cold.Ny=pipe_cold.Eta/pipe_cold.density; %кинематическая вязкость

%pipe_cold.k=(2*pipe_cold.lambda)/(pipe_cold.douter*log(pipe_cold.douter/pipe_cold.din)); % коэффициент теплопередачи

end