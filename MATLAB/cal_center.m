function center = cal_center(A, B, C)
% 根据三个空间点，计算出其圆心及半径

center = 0;
% 数据检查
% 检查数据输入格式是否正确
if size(A,2)~=3 || size(B,2)~=3 || size(C,2)~=3
    fprintf('输入点维度不一致\n');return;
end
n = size(A,1);
if size(B,1)~=n || size(C,1)~=n
    fprintf('输入点维度不一致\n');return;
end

% 计算A到B的单位向量和A到C的单位向量
% 检查点是否相同
ab = B - A;
ac = C - A;
if find(norm(ab)==0) | find(norm(ac)==0) %#ok<OR2>
    fprintf('输入点不能一样\n');return;
end

ab_1 = ab/norm(ab);
ac_1 = ac/norm(ac);

% 计算圆平面上的单位法向量
% 检查三点是否共线
n = cross(ab_1,ac_1);
 if all(n==0)
    fprintf('三个点共线\n');return;
 end
if find(sum(abs(n),2)<1e-5)
    fprintf('三点过于趋近直线\n');return;
end

% 计算新坐标系UVW轴
u = ab_1;
w = cross(ac,ab)/norm(cross(ac,ab));
v = cross(w,u);

% 计算投影
bx = dot(ab,u);
cx = dot(ac,u);
cy = dot(ac,v);

% 计算圆心
h = ((cx - bx/2)^2 + cy^2 -(bx/2)^2)/(2*cy);
center = zeros(1,3);
center(1,:) = A(1,:) + bx/2.*u(1,:) + h.*v(1,:);

end
