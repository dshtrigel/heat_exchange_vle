function NewtonCounterFlow()
clc;
close all;

% Входные параметры
Thot_in = 173;
Tcold_out = 166.1;

% Считаем вспомогательную задачу. 
% Для заданной выходной температуры холодного потока ищем соответствующую
% входную. Мы тем самым генерируем исходную информацию, то, что при
% боевом расчете нам изначально известно
[~, Tcold_in] = CalcEulerCounterFlowCont(Thot_in, Tcold_out);

% создаем слой абстракции в виде уравнения f(x) = 0
f = @(x) CalcEulerCounterFlow_Tcold_in(Thot_in, x) - Tcold_in;

[Tcold_out_calc, n] = NewtonMethod(f, 150); % должно совпасть с Tcold_out
Tcold_out_calc
Tcold_out
end

% Вспомогательная функция, вытаскивает только Tcold_in из решения
% вспомогательной задачи
function Tcold_in = CalcEulerCounterFlow_Tcold_in(Thot_in, Tcold_out)
[~, Tcold_in] = CalcEulerCounterFlowCont(Thot_in, Tcold_out);

end

