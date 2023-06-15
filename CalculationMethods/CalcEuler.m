% Метод Эйлера
function [Tx_hot, Tx_cold]=CalcEuler(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold, type)

S_hot= pi * pipe_hot.din^2 / 4;
v_hot=Ghot/(pipe_hot.density*S_hot);
Eps_hot=0.015/pipe_hot.din;
Re_hot=v_hot*pipe_hot.din/pipe_hot.Ny;
lambdatr_hot=CalcLambdaAltshul(Re_hot, Eps_hot);

S_cold = pi * pipe_cold.din^2 / 4;
v_cold=Gcold/(pipe_cold.density*S_cold);
Eps_cold=0.015/pipe_cold.din;
Re_cold=v_cold*pipe_cold.din/pipe_cold.Ny;
lambdatr_cold=CalcLambdaAltshul(Re_cold, Eps_cold);

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);

Tx_hot(1)=Thot_in;
Tx_cold(1)=Tcold_in;
n=length(x);
for i=1:n-1
    t_hot = Tx_hot(i);
    t_cold = Tx_cold(i);
   
    q_hot = -k * (t_hot-t_cold);
    q_cold= k * (t_hot-t_cold);

    P=pi*(pipe_hot.din+pipe_cold.din)/2;

    if type == "fric"
        f_hot=((q_hot)*P/(Ghot*pipe_hot.c))+((lambdatr_hot*v_hot^2)/(2*pipe_hot.din*pipe_hot.c));
        f_cold=((q_cold)*P/(Gcold*pipe_cold.c))+((lambdatr_cold*v_cold^2)/(2*pipe_cold.din*pipe_cold.c));
    elseif type == "without_fric"
        f_hot=((q_hot)*P/(Ghot*pipe_hot.c));
        f_cold=((q_cold)*P/(Gcold*pipe_cold.c));
    else
        exception = MException('foo:BadInput', 'Not enough inputs');
        throw(exception);
    end
    dx = x(i+1) - x(i);
    Tx_hot(i+1)=t_hot+dx*f_hot;
    Tx_cold(i+1)=t_cold+dx*f_cold;
end 
end