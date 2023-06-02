function [x, i] = NewtonMethod(func, x0)
    % функция метода Ньютона
    % func - функция, которую нужно решить
    % df - производная функции
    % TOutCold0 - начальное приближение
    % TInHot - температура горячего потока на входе
    % TInColdSet - заданная температура холодного потока на входе
    eps = 0.001; % допустимая погрешность
    max_iter = 100; % максимальное количество итераций

    % Задаем начальное значение
    x = x0;

    for i = 1:max_iter
        % Вычисляем значение функции и ее производной в точке TOutCold
        f = func(x);
        dfx = df(func, x);

        % Обновляем TOutCold с использованием формулы Ньютона
        dx = -f/dfx;
        x = x + dx;% новое значение переменной

        % Проверяем, достигли ли мы точности
        if abs(dx) < eps
            break;
        end
    end
end

function dfx = df(f, x)
    h = 0.01;  % small step
    f_plus = f(x + h);
    f_minus = f(x - h);
    dfx = (f_plus - f_minus) / (2*h);
end
