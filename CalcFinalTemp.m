function CalcFinalTemp()
clc;
close all;

% Входные параметры
pipe_hot = GetPipe_hot();
pipe_cold = GetPipe_cold();


x = 0:0.01:pipe_hot.L;

% Расчет конечных температур методом Эйлера с трением
[Tx_hot_fric, Tx_cold_fric]=CalcEuler(173, 153, 50, 35, x, pipe_hot, pipe_cold, "fric");
% Расчет конечных температур методом Эйлера без трения
[Tx_hot_NoFric, Tx_cold_NoFric]=CalcEuler(173, 153, 50, 35, x, pipe_hot, pipe_cold, "without_fric");
% Расчет конечных температур аналитическим методом (без трения)
[Thot_out, Tcold_out]= CalcAnalytFinalTemp(173, 153, 50, 35, x, pipe_hot, pipe_cold);

% Вывод графиков
figure
plot(x,Tx_hot_fric, 'LineWidth', 1)
hold on
plot(x,Tx_cold_fric, 'LineWidth', 1)
hold on
plot(x,Thot_out, 'LineWidth', 3)
hold on
plot(x,Tcold_out, 'LineWidth', 3)
hold on
plot (x,Tx_hot_NoFric, 'LineWidth', 1); 
hold on
plot (x,Tx_cold_NoFric, 'LineWidth', 1);
hold on
legend('T_x hot fric', 'T_x cold fric', 'T hot out', 'T cold out', 'T_x hot NoFric', 'T_x cold NoFric', 'Location','southeast')
hold off

end