function center = cal_center(A, B, C)
% ���������ռ�㣬�������Բ�ļ��뾶

center = 0;
% ���ݼ��
% ������������ʽ�Ƿ���ȷ
if size(A,2)~=3 || size(B,2)~=3 || size(C,2)~=3
    fprintf('�����ά�Ȳ�һ��\n');return;
end
n = size(A,1);
if size(B,1)~=n || size(C,1)~=n
    fprintf('�����ά�Ȳ�һ��\n');return;
end

% ����A��B�ĵ�λ������A��C�ĵ�λ����
% �����Ƿ���ͬ
ab = B - A;
ac = C - A;
if find(norm(ab)==0) | find(norm(ac)==0) %#ok<OR2>
    fprintf('����㲻��һ��\n');return;
end

ab_1 = ab/norm(ab);
ac_1 = ac/norm(ac);

% ����Բƽ���ϵĵ�λ������
% ��������Ƿ���
n = cross(ab_1,ac_1);
 if all(n==0)
    fprintf('�����㹲��\n');return;
 end
if find(sum(abs(n),2)<1e-5)
    fprintf('�����������ֱ��\n');return;
end

% ����������ϵUVW��
u = ab_1;
w = cross(ac,ab)/norm(cross(ac,ab));
v = cross(w,u);

% ����ͶӰ
bx = dot(ab,u);
cx = dot(ac,u);
cy = dot(ac,v);

% ����Բ��
h = ((cx - bx/2)^2 + cy^2 -(bx/2)^2)/(2*cy);
center = zeros(1,3);
center(1,:) = A(1,:) + bx/2.*u(1,:) + h.*v(1,:);

end
