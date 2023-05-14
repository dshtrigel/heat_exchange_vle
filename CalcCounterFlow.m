function CalcCounterFlow()
clc;
close all;

% Входные параметры
pipe_hot = GetPipe_hot();
pipe_cold = GetPipe_cold();
Thot_in = 173;
Tcold_out = 166.1;

x = 0:0.01:pipe_hot.L;


[Tx_hot, Tx_cold]=CalcEulerCounterFlow(Thot_in, Tcold_out, 50, 35, x, pipe_hot, pipe_cold);

% TInHot = Tx_hot(1)
%TOutCold = Tx_cold(1)
% 
% TOutHot = Tx_hot(length(x))
% TInCold = Tx_cold(length(x))


[Thot_out, Tcold_in] = CalcEulerCounterFlowCont(Thot_in, Tcold_out);
Thot_out
Tcold_in
Tcold_out_calc = Newton(Thot_in, Tcold_in, 150)

plot (x,Tx_hot, x,Tx_cold)
legend ('T_x hot', 'T_x cold')
end



