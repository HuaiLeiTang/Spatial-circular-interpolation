%��ʼ��ʱ
tic
clear

%��ʼ�㡢��ֹ��

A = [0.3,0.15,0.32];
C = [0.3,0.15,0.89];

%����dz
deta_d = cal_d(A,C);

%�����м��
B = cal_pointB(A,C,deta_d);

%ˮƽ��������
%if d(3) == 0
   % p2 = CalP2_0(p3,p1);
%��ֱ��xoyƽ��
%elseif d(2) == 0 && d(1) == 0
    %p2 = CalP2_1(p3,p1);
%else
    %p2 = CalP2_2(p3,p1);
%end

%�岹����
stepL = cal_stepL(A,C);

center = cal_center(A, B, C);
rad = cal_rad(A,center);

%�岹
points = cal_Inter(A,B,C,center,rad,stepL);

num = size(points,2);
for i=1:num
    plot3(points(1,i),points(2,i),points(3,i),'-*r');
    hold on
end

axis ([-1 1 -1 1 0 1])
grid on

%������ʱ
toc