function ClacTubeHeatExchanger()
clc;
close all;
pipe_hot = GetPipe_hot();
pipe_cold = GetPipe_cold();


x = 0:0.01:pipe_hot.L;

%[Tx_hot, Tx_cold]=CalcEuler(173, 153, 50, 35, x, pipe_hot, pipe_cold);
[Thot_out, Tcold_out]= FinalTemperatures(173, 153, 50, 35, x, pipe_hot, pipe_cold);
[Tx_hot2, Tx_cold2]=CalcEuler3(173, 153, 50, 35, x, pipe_hot, pipe_cold);
% [Thot_out2, Tcold_out2]= FinalTemperatures2(173, 153, 50, 35, x, pipe_hot, pipe_cold); 


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
legend('T hot out Analyt', 'T cold out Analyt', 'T_x hot EulerNoFric', 'T_x cold EulerNoFric')
hold off

end



% Метод Эйлера
function [Tx_hot, Tx_cold]=CalcEuler(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

S_hot= pi * pipe_hot.din^2 / 4;
v_hot=Ghot/(pipe_hot.density*S_hot);
Eps_hot=0.015/pipe_hot.din;
Re_hot=v_hot*pipe_hot.din/pipe_hot.Ny;
lambdatr_hot=CalcLambdaAltshul(Re_hot, Eps_hot);

S_cold = pi * pipe_cold.din^2 / 4;
v_cold=Gcold/(pipe_cold.density*S_cold);
Eps_cold=0.015/pipe_cold.din;
Re_cold=v_cold*pipe_cold.din/pipe_cold.Ny;
lambdatr_cold=CalcLambdaAltshul(Re_cold, Eps_cold);

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + 1/pipe_cold.alpha);

Tx_hot(1)=Thot_in;
Tx_cold(1)=Tcold_in;
n=length(x);
for i=1:n-1
    t_hot = Tx_hot(i);
    t_cold = Tx_cold(i);
   
    q_hot = -k * (t_hot-t_cold);
    q_cold= k * (t_hot-t_cold);

    P=pi*pipe_cold.douter;

    f_hot=((q_hot)*P/(Ghot*pipe_hot.c))+((lambdatr_hot*v_hot^2)/(2*pipe_hot.din*pipe_hot.c));
    f_cold=((q_cold)*P/(Gcold*pipe_cold.c))+((lambdatr_cold*v_cold^2)/(2*pipe_cold.din*pipe_cold.c));

    dx = x(i+1) - x(i);
    Tx_hot(i+1)=t_hot+dx*f_hot;
    Tx_cold(i+1)=t_cold+dx*f_cold;
end 
end

% метод Эйлера без трения
function [Tx_hot2, Tx_cold2]=CalcEuler3(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

S_hot= pi * pipe_hot.din^2 / 4;
v_hot=Ghot/(pipe_hot.density*S_hot);

S_cold = pi * pipe_cold.din^2 / 4;
v_cold=Gcold/(pipe_cold.density*S_cold);

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + 1/pipe_cold.alpha);

Tx_hot2(1)=Thot_in;
Tx_cold2(1)=Tcold_in;
n=length(x);
for i=1:n-1
    t_hot = Tx_hot2(i);
    t_cold = Tx_cold2(i);
    q_hot = -k * (t_hot-t_cold);
    q_cold= k * (t_hot-t_cold);

    P=pi*pipe_cold.douter;

    f_hot=(q_hot*P/(Ghot*pipe_hot.c));
    f_cold=(q_cold*P/(Gcold*pipe_cold.c));

%     f_hot=(-(4*q_hot)/(pipe_hot.din*pipe_hot.density*pipe_hot.c*v_hot));
%     f_cold=(-(4*q_cold)/(pipe_cold.din*pipe_cold.density*pipe_cold.c*v_cold));
    
    dx = x(i+1) - x(i);
    Tx_hot2(i+1)=t_hot+dx*f_hot;
    Tx_cold2(i+1)=t_cold+dx*f_cold;
end 
end


% вывод по Исаченко
function [Thot_out, Tcold_out]= FinalTemperatures(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + 1/pipe_cold.alpha);

Chot=pipe_hot.c*Ghot;
Ccold=pipe_cold.c*Gcold;
m = 1/Chot+1/Ccold;

F=pi*x*pipe_cold.douter;
Thot_out=Thot_in-(Thot_in-Tcold_in)*(1-exp(-k*F*m))/(1+Chot/Ccold);
Tcold_out=Tcold_in+(Thot_in-Tcold_in)*(1-exp(-k*F*m))/(1+Ccold/Chot);

end
