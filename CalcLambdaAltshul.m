% Определение коэффициента гидравлического трения 
function lambdatr = CalcLambdaAltshul(Re, Eps)
lambdatr=0.11*(Eps+68/Re)^(1/4);
end