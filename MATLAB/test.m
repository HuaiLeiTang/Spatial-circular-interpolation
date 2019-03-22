%开始计时
tic
clear

%起始点、终止点

A = [0.3,0.15,0.32];
C = [0.3,0.15,0.89];

%计算dz
deta_d = cal_d(A,C);

%计算中间点
B = cal_pointB(A,C,deta_d);

%水平的两个点
%if d(3) == 0
   % p2 = CalP2_0(p3,p1);
%垂直于xoy平面
%elseif d(2) == 0 && d(1) == 0
    %p2 = CalP2_1(p3,p1);
%else
    %p2 = CalP2_2(p3,p1);
%end

%插补次数
stepL = cal_stepL(A,C);

center = cal_center(A, B, C);
rad = cal_rad(A,center);

%插补
points = cal_Inter(A,B,C,center,rad,stepL);

num = size(points,2);
for i=1:num
    plot3(points(1,i),points(2,i),points(3,i),'-*r');
    hold on
end

axis ([-1 1 -1 1 0 1])
grid on

%结束计时
toc