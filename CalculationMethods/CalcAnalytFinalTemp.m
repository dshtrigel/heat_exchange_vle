% Аналитические способы определения конечных температур
function [Thot_out, Tcold_out]= CalcAnalytFinalTemp(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)

k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);

Chot=pipe_hot.c*Ghot;
Ccold=pipe_cold.c*Gcold;
m = 1/Chot+1/Ccold;

F=pi*x*(pipe_hot.din+pipe_cold.din)/2;
Thot_out=Thot_in-(Thot_in-Tcold_in)*(1-exp(-k*F*m))/(1+Chot/Ccold);
Tcold_out=Tcold_in+(Thot_in-Tcold_in)*(1-exp(-k*F*m))/(1+Ccold/Chot);

end

% % По формуле Белоконя
% function [Thot_out2, Tcold_out2]= FinalTemperatures2(Thot_in, Tcold_in, Ghot, Gcold, x, pipe_hot, pipe_cold)
% 
% k = 1 /(1/pipe_hot.alpha + pipe_hot.wall/pipe_hot.lambda + pipe_cold.wall/pipe_cold.lambda + 1/pipe_cold.alpha);
% Chot=pipe_hot.c*Ghot;
% Ccold=pipe_cold.c*Gcold;
% C = ((Chot*Ccold)/(Chot+Ccold));
% n=length(x);
% for i=1:n
%     F = pi*x(i)*(pipe_hot.din+pipe_cold.din)/2; 
%     if x(i)==0
%        Q=0;
%     else
%        Q=2*(Thot_in-Tcold_in)/(1/Chot + 1/Ccold + (1/C)*((exp(k*F/C)+1)/(exp(k*F/C)-1)));  
%     end
%     Thot_out2(i)=Thot_in-Q/Chot;
%     Tcold_out2(i)=Tcold_in+Q/Ccold;
% end
% end   
% 
% % n=length(x);
% % 
% % Q = zeros(size(x));
% % Thot_out2 = zeros(size(x));
% % Tcold_out2 = zeros(size(x));
% % 
% % % Заполняется 1 элемент необходимых массивов
% % Q(:, 1) = 0;
% % Thot_out2(:, 1) = Thot_in;
% % Tcold_out2(:, 1) = Tcold_in;
% % 
% % % Заполенение всех остальных 2:n элементов
% % for i=2:n
% %     F = pi*x(:, i)*(pipe_hot.douter+pipe_cold.douter)/2;
% %     Q(:, i)=2*(Thot_in-Tcold_in)/(1/Chot + 1/Ccold) + (1/C)*((exp(k*F/C)+1)/(exp(k*F/C)-1)));  
% %     Thot_out2(:, i)=Thot_in-Q(:, i)/Chot;
% %     Tcold_out2(:, i)=Tcold_in+Q(:, i)/Ccold;
% %    
% % end
% %end