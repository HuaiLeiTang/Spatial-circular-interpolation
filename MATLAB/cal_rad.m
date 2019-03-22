function rad = cal_rad(A,center)
% 根据空间点、圆心，计算出其半径

% 半径
rad = sqrt((center(1,1)-A(1,1)).^2+(center(1,2)-A(1,2)).^2+(center(1,3)-A(1,3)).^2);
end
