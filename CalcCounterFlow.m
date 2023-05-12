function CalcCounterFlow()
clc;
close all;

% Входные параметры
pipe_hot = GetPipe_hot();
pipe_cold = GetPipe_cold();

x = 0:0.01:pipe_hot.L;

[Tx_hot, Tx_cold]=CalcEulerCounterFlow(173, 166.1, 50, 35, x, pipe_hot, pipe_cold);

plot (x,Tx_hot, x,Tx_cold)
legend ('T_x hot', 'T_x cold')
end



