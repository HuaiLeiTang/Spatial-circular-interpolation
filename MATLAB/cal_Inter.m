function points = cal_Inter(A,B,C,center,rad,stepL)

% 建立圆弧坐标系
% 计算转换矩阵
L = (B(2)-A(2))*(C(3)-B(3))-(B(3)-A(3))*(C(2)-B(2));
M = (B(3)-A(3))*(C(1)-B(1))-(B(1)-A(1))*(C(3)-B(3));
N = (B(1)-A(1))*(C(2)-B(2))-(B(2)-A(2))*(C(1)-B(1));
K = sqrt(L^2+M^2+N^2);
a = [L M N]/K;
n = (A -center)/rad;
o = cross(a,n);
T = [n' o' a' center'; 0 0 0 1];

% 求转换后的点
q1 = inv(T)*[A 1]';
q2 = inv(T)*[B 1]';
q3 = inv(T)*[C 1]';

% 计算角度
if q3(2)<0
    theta13 = atan2(q3(2),q3(1)) + 2*pi;
else
    theta13 = atan2(q3(2),q3(1));
end

if q2(2)<0
    theta12 = atan2(q2(2),q2(1)) + 2*pi;
else
    theta12 = atan2(q2(2),q2(1));
end

% 轨迹插补
count =1;

for step = 0:stepL:theta13
    points(:,count) = T*[rad*cos(step) rad*sin(step) 0 1]';
    count = count+1;
    residual = theta13 - step;
    if residual <= stepL
        break;
    end
end

%p(:,count) = T*[rad*cos(residual) rad*sin(residual) 0 1]';
points(:,count) = [C(1) C(2) C(3) 1]';

end
