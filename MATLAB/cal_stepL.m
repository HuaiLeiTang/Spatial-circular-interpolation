function stepL = cal_stepL(A,C)
% 根据起始、终止点计算出步长

param = 0.5;

%p1、p3的距离
distance = sqrt((C(1)-A(1)).^2+(C(2)-A(2)).^2+(C(3)-A(3)).^2);

if distance <= 2*param
    stepL = 0.04;
elseif distance <= 3*param
    stepL = 0.02;
elseif distance <= 4*param
    stepL = 0.012;
else
    stepL = 0.008;
end

end
    