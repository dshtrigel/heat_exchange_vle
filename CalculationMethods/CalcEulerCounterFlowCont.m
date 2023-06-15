% метод Эйлера без трения для противотока
function [Thot_out, Tcold_in]=CalcEulerCounterFlowCont(Thot_in, Tcold_out)

Ghot = 50;
Gcold = 35;

pipe_hot = GetPipe_hot();
pipe_cold = GetPipe_cold();

x = 0:0.01:pipe_hot.L;

S_hot= pi * pipe_hot.din^2 / 4;
v_hot=Ghot/(pipe_hot.density*S_hot);

S_cold = pi * pipe_cold.din^2 / 4;
v_cold=Gcold/(pipe_cold.density*S_cold);

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);

Tx_hot(1)=Thot_in;
Tx_cold(1)=Tcold_out;
n=length(x);
for i=1:n-1
    t_hot = Tx_hot(i);
    t_cold = Tx_cold(i);
    q_hot = -k * (t_hot-t_cold);
    q_cold= k * (t_hot-t_cold);

    P=pi*(pipe_hot.din+pipe_cold.din)/2;

    f_hot=(q_hot*P/(Ghot*pipe_hot.c));
    f_cold=(-q_cold*P/(Gcold*pipe_cold.c));
    
    dx = x(i+1) - x(i);
        Tx_hot(i+1)=t_hot+dx*f_hot;
        Tx_cold(i+1)=t_cold+dx*f_cold;
end 
Thot_out = Tx_hot(n);
Tcold_in =  Tx_cold(n);
end