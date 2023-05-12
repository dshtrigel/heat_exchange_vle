% Формирование параметров горячей трубы
function pipe_hot = GetPipe_hot()
pipe_hot.L=50; % длина

pipe_hot.din=0.1; % внутренний диаметр
pipe_hot.douter=0.108; % внешний диаметр
pipe_hot.wall=0.004; % толщина стенки

pipe_hot.lambda=209; % коэффициент теплопроводности

pipe_hot.alpha=52250;% коэффициент теплоотдачи
pipe_hot.c=5477; % теплоемкость
pipe_hot.density = 419; % плотность СПГ
pipe_hot.Eta=0.0000049;% динамическая вязкость
pipe_hot.Ny=pipe_hot.Eta/pipe_hot.density; %кинематическая вязкость

pipe_hot.k=(2*pipe_hot.lambda)/(pipe_hot.douter*log(pipe_hot.douter/pipe_hot.din)); % коэффициент теплопередачи

end