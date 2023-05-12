% Формирование параметров холодной трубы
function pipe_cold = GetPipe_cold()
pipe_cold.L=50; % длина

pipe_cold.din=0.08; % внутренний диаметр
pipe_cold.douter=0.088; % внешний диаметр
pipe_cold.wall=0.004; % толщина стенки

pipe_cold.lambda=209; % коэффициент теплопроводности

pipe_cold.alpha=52250; % коэффициент теплоотдачи
pipe_cold.c=3950; % теплоемкость
pipe_cold.density =1349; % плотность азота
pipe_cold.Eta=0.0000129;% динамическая вязкость
pipe_cold.Ny=pipe_cold.Eta/pipe_cold.density; %кинематическая вязкость

pipe_cold.k=(2*pipe_cold.lambda)/(pipe_cold.douter*log(pipe_cold.douter/pipe_cold.din)); % коэффициент теплопередачи

end