function rad = cal_rad(A,center)
% ���ݿռ�㡢Բ�ģ��������뾶

% �뾶
rad = sqrt((center(1,1)-A(1,1)).^2+(center(1,2)-A(1,2)).^2+(center(1,3)-A(1,3)).^2);
end
