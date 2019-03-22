function p2 = CalP2_2(p1,p3)
% 根据起始、终止点计算出中间点

p2(1)=(p1(1)+p3(1))/2;
p2(2)=(p1(2)+p3(2))/2;

%p1、p3的距离
%d1 = sqrt((p3(1)-p1(1)).^2+(p3(2)-p1(2)).^2+(p3(3)-p1(3)).^2);
dz = abs(p3(3)-p1(3));

if dz <= 0.3
    h = dz/3;
elseif dz <= 0.5
    h = dz/4;
else
    h = dz/6;
end

p2(3)=(p1(3)+p3(3))/2 + h;

end