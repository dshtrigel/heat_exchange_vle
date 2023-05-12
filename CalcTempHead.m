function CalcTempHead()
clc;
close all;

% Входные параметры
pipe_hot = GetPipe_hot();
pipe_cold = GetPipe_cold();


x = 0:0.1:pipe_hot.L;

dT_Euler = CalcEulerTempHead(173, 153, 50, 35, x, pipe_hot, pipe_cold);
dT_Analyt = CalcAnalytTempHead(173, 153, 50, 35, x, pipe_hot, pipe_cold);

figure
plot(x,dT_Euler, 'LineWidth', 1.5)
hold on
plot(x, dT_Analyt, 'LineWidth', 1)
hold on
legend('dT Euler', 'dT Analyt')
hold off
end





