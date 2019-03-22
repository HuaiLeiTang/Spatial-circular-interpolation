function B = cal_pointB(A,C,deta_d)
% ������ʼ����ֹ�������м��

param=0.5;

%p1��p3�ľ���
distance = sqrt((C(1)-A(1)).^2+(C(2)-A(2)).^2+(C(3)-A(3)).^2);

if deta_d(3)==0   %ˮƽ��������
    B(1)=(A(1)+C(1))/2;
    B(2)=(A(2)+C(2))/2;
    
    if distance <= param
        h = distance/12;
    elseif distance <= 2*param
        h = distance/20;
    elseif distance <= 3*param
        h = distance/25;
    else
        h = distance/30;
    end
    
    B(3)=(A(3)+C(3))/2 + h;
    
elseif deta_d(2) == 0 && deta_d(1) == 0       %AC��ֱ��xoyƽ��
    B(1)=(A(1)+C(1))/2;
    B(3)=(A(3)+C(3))/2;
    
    if distance <= param
        h = distance/12;
    elseif distance <= 2*param
        h = distance/20;
    elseif distance <= 3*param
        h = distance/25;
    else
        h = distance/30;
    end
    
    B(2)=(A(2)+C(2))/2 + h;
    
else    %��ͨ���
    B(1)=(A(1)+C(1))/2;
    B(2)=(A(2)+C(2))/2;
    
    if deta_d(3) <= 0.3
        h = deta_d(3)/3;
    elseif deta_d(3) <= 0.5
        h = deta_d(3)/4;
    else
        h = deta_d(3)/6;
    end

    B(3)=(A(3)+C(3))/2 + h;
    
end