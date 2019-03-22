function deta_d = cal_d(A,C)
%º∆À„dz

dx = abs(C(1)-A(1));
dy = abs(C(2)-A(2));
dz = abs(C(3)-A(3));

deta_d = [dx,dy,dz];

end