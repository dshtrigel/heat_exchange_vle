function [output] = Newton(input, set_value, init)
    h = 0.001; % приращение независимого аргумента при вычислении производной
   hinv2 = 1/(2*h); % множитель для вычисления производной
    eps_max = 0.01; % максимально допусимая ошибка по x
    dx = 1d+12; % стартовое значение ошибки
    x = init;% начальное приближение
    niter = 0; % число итераций устанавливается на ноль
    while abs(dx) > eps_max % пока ошибка превышает допустимую,
        niter = niter + 1; % увеличить счётчик числа итераций на 1
        y = CalcEulerCounterFlowCont(input, x); % значение функции в текущей точке
        dy = f1(x,h,hinv2,input, set_value); % значение производной в текущей точке
        dx = y/dy; % смещение по x для перехода ближе к корню
        x = x - dx; % новое значение переменной
%        disp(x)
    end % конец цикла
     disp('выполнено итераций:'); 
     disp(niter); % вывести число итераций
    x = round_2(x); % округление корня до сотых
%     disp('найденное значение выходной температуры холодного потока:');
%     disp(x); % показать значение корня
    output = x;
end

function dy1=f1(x1,hloc,hinvloc, input, set_value) % подпрограмма вычисления производной
    dy1 = hinvloc*(func(input, x1+hloc, set_value)-func(input, x1-hloc, set_value)); % значение производной в точке
end % конец подпрограммы вычисления производной

function u=round_2(u) % подпрограмма округления до сотых
    u = 0.01*round(100*u); % умножить на 100, округлить до целых и разделить на 100 
end % конец подпрограммы

function [delta_y]=func(Thot_in, x, set_value)
[Thot_out, Tcold_in] = CalcEulerCounterFlowCont(Thot_in, x);
delta_y =  Tcold_in - set_value; 
end % Конец основной программы