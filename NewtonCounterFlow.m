function NewtonCounterFlow()
clc;
close all;

function TOutCold = newtonMethod(func, df, TOutCold0, TInHot, TInColdSet)
    % функция метода Ньютона
    % func - функция, которую нужно решить
    % df - производная функции
    % TOutCold0 - начальное приближение
    % TInHot - температура горячего потока на входе
    % TInColdSet - заданная температура холодного потока на входе
    tol = 0.001; % допустимая погрешность
    max_iter = 100; % максимальное количество итераций

    % Задаем начальное значение
    TOutCold = TOutCold0;

    for i = 1:max_iter
        % Вычисляем значение функции и ее производной в точке TOutCold
        [~,TInColdCalc] = func(TInHot, TOutCold);
        f = TInColdCalc - TInColdSet;
        dfx = df(TInHot, TOutCold);

        % Обновляем TOutCold с использованием формулы Ньютона
        dT = -f/dfx;
        TOutCold = TOutCold + dT;

        % Проверяем, достигли ли мы точности
        if abs(dT) < tol
            break;
        end
    end
end

function dfx = df(TInHot, TOutCold)
    h = 0.01;  % small step
    [~ , TInColdCalc1] = CalcEulerCounterFlowCont(TInHot, TOutCold + h);
    [~ , TInColdCalc2] = CalcEulerCounterFlowCont(TInHot, TOutCold - h);
    dfx = (TInColdCalc1 - TInColdCalc2) / (2*h);
end

% Входные параметры
Thot_in = 173;
Tcold_out = 166.1;

[Thot_out, Tcold_in] = CalcEulerCounterFlowCont(Thot_in, Tcold_out);

TInColdSet = Tcold_in; % - заданная температура холодного потока на входе

TOutCold = newtonMethod(@CalcEulerCounterFlowCont, @df, 150, Thot_in, TInColdSet)

end