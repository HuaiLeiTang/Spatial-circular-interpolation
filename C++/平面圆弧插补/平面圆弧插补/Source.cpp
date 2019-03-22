#include "time.h"
#include "math.h"
#include "stdio.h"
#include "iostream"
#include "Eigen/Dense"
#include "Eigen/Geometry"
//#include "circular_interpolation.h"

#define PI acos(-1)

using namespace std;
using namespace Eigen;

// circular interpolation
class circle
{
public:
	//RowVector3d deta_d;
	double distance;
	double param = 0.5;
	double cal_stepL(RowVector3d a, RowVector3d b);
	RowVector3d cal_dz(RowVector3d a, RowVector3d b);
	RowVector3d cal_pointB(RowVector3d a, RowVector3d b, RowVector3d deta_d);
	RowVector3d cal_center(RowVector3d a, RowVector3d b, RowVector3d c, double *rad);
	MatrixXd  cal_Inter(RowVector3d a, RowVector3d b, RowVector3d c, RowVector3d center, double rad, double stepL);
};

RowVector3d circle::cal_dz(RowVector3d a, RowVector3d b)
{
	RowVector3d deta_d;
	double dx, dy, dz;

	dx = abs(a(0) - b(0));
	dy = abs(a(1) - b(1));
	dz = abs(a(2) - b(2));

	deta_d(0) = dx;
	deta_d(1) = dy;
	deta_d(2) = dz;

	return deta_d;
}

RowVector3d circle::cal_pointB(RowVector3d a, RowVector3d b, RowVector3d deta_d)
{
	double h;
	RowVector3d pointB;

	distance = sqrt(pow(deta_d(0), 2) + pow(deta_d(1), 2) + pow(deta_d(2), 2));

	if (deta_d(2) == 0)    //A and B are horizontal
	{
		pointB(0) = (a(0) + b(0)) / 2;
		pointB(1) = (a(1) + b(1)) / 2;

		if (distance <= param)
			h = distance / 12;
		else if (distance <= 2 * param)
			h = distance / 20;
		else if (distance <= 3 * param)
			h = distance / 25;
		else
			h = distance / 30;

		pointB(2) = (a(2) + b(2)) / 2 + h;
	}

	else if ((deta_d(1) == 0) && (deta_d(0) == 0))    //Vertical to the plane xoy
	{
		pointB(0) = (a(0) + b(0)) / 2;
		pointB(2) = (a(2) + b(2)) / 2;

		if (distance <= param)
			h = distance / 12;
		else if (distance <= 2 * param)
			h = distance / 20;
		else if (distance <= 3 * param)
			h = distance / 25;
		else
			h = distance / 30;

		pointB(1) = (a(1) + b(1)) / 2 + h;
	}

	else   //Generally
	{
		pointB(0) = (a(0) + b(0)) / 2;
		pointB(1) = (a(1) + b(1)) / 2;

		if (deta_d(2) <= 0.3)
			h = deta_d(2) / 3;
		else if (deta_d(2) <= 0.5)
			h = deta_d(2) / 4;
		else
			h = deta_d(2) / 6;

		pointB(2) = (a(2) + b(2)) / 2 + h;
	}

	return pointB;
}

double circle::cal_stepL(RowVector3d a, RowVector3d b)
{
	double stepL;

	if (distance <= 2 * param)
		stepL = 0.04;
	else if (distance <= 3 * param)
		stepL = 0.02;
	else if (distance <= 4 * param)
		stepL = 0.012;
	else
		stepL = 0.008;

	return stepL;
}

RowVector3d circle::cal_center(RowVector3d a, RowVector3d b, RowVector3d c, double *rad)
{
	RowVector3d center, ab, ac, n;
	RowVector3d u, v, w;
	double bx, cx, cy, center_uvw;

	//vector ab
	ab = b - a;
	ac = c - a;
	n = ac.cross(ab);       //叉乘：cross(ac,ab) ,requires #include "Eigen/Geometry"

	//Calculate the new coordinate system UVW axis
	u = ab / ab.norm();
	w = n / n.norm();
	v = w.cross(u);

	//a在UVW坐标系中的位置为(0,0,0);b为(bx,0,0);c为(cx,cy,0)
	//Calculate the projection on the UVW axis
	bx = ab.dot(u);    //点积：dot(ab,u)
	cx = ac.dot(u);
	cy = ac.dot(v);
	
	//设在UVW坐标系中，圆心坐标为(bx/2,center_uvw,0)
	center_uvw = (pow((cx - bx / 2), 2) + pow(cy, 2) - pow(bx / 2, 2)) / (2 * cy);
	center = a + (bx / 2)*u + center_uvw*v;

	*rad = sqrt(pow((center(0) - a(0)), 2) + pow((center(1) - a(1)), 2) + pow((center(2) - a(2)), 2));
	return center;
}

MatrixXd circle::cal_Inter(RowVector3d a, RowVector3d b, RowVector3d c, RowVector3d center, double rad, double stepL)
{
	Matrix4d T;
	RowVector3d O, P, Q;
	Vector4d a_1, b_1, c_1, a_uvw, b_uvw, c_uvw, step;
	double l, m, n, k;
	double theta_ac, theta_ab,residual;
	double alpha = 0;
	int count, i;

	//计算转换矩阵
	l = (b(1) - a(1))*(c(2) - b(2)) - (b(2) - a(2))*(c(1) - b(1));
	m = (b(2) - a(2))*(c(0) - b(0)) - (b(0) - a(0))*(c(2) - b(2));
	n = (b(0) - a(0))*(c(1) - b(1)) - (b(1) - a(1))*(c(0) - b(0));
	k = sqrt(pow(l, 2) + pow(m, 2) + pow(n, 2));

	O(0) = l / k;            //a
	O(1) = m / k;            
	O(2) = n / k;            
	P = (a - center) / rad;  //n
	Q = O.cross(P);          //o

	T.block(0, 0, 3, 1) = P.transpose();
	T.block(0, 1, 3, 1) = Q.transpose();
	T.block(0, 2, 3, 1) = O.transpose();
	T.block(0, 3, 3, 1) = center.transpose();
	T.block(3, 0, 1, 4) << 0, 0, 0, 1;

	//求A、B、C点在新坐标系UVW下对应的坐标值
	a_1 << a(0), a(1), a(2), 1;
	b_1 << b(0), b(1), b(2), 1;
	c_1 << c(0), c(1), c(2), 1;

	a_uvw = T.inverse()*a_1;
	b_uvw = T.inverse()*b_1;
	c_uvw = T.inverse()*c_1;

	//计算角度
	if (c_uvw(1) < 0)
		theta_ac = atan2(c_uvw(1), c_uvw(0)) + 2 * PI;
	else
		theta_ac = atan2(c_uvw(1), c_uvw(0));
	
	if (b_uvw(1) < 0)
		theta_ab = atan2(b_uvw(1), b_uvw(0)) + 2 * PI;
	else
		theta_ab = atan2(b_uvw(1), b_uvw(0));

	//轨迹插补
	count = theta_ac / stepL;
	MatrixXd points(4, count+2);
	points.col(0) << a(0), a(1), a(2), 1;
	for (i = 1; i < count+1; i++)
	{
		alpha = alpha + stepL;
		step << rad*cos(alpha), rad*sin(alpha), 0, 1;
		points.col(i) = T*step;
	}
	//residual = theta_ac - count*stepL;   //剩余的角度
	 points.col(count+1) << c(0), c(1), c(2), 1;

	return points;
}


int main()
{
	//Calulate time
	clock_t start, finish;
	double duration;
	start = clock();

	circle inter;
	//RowVector3d A, B, C;
	RowVector3d B;
	RowVector3d d_xyz, center;
	double stepL, rad;
	MatrixXd points;
	int i, j;

	/*//Input A、C
	cout << "Input start point (A):" << endl;
	cin >> A(0) >> A(1) >> A(2);

	cout << "Input end point (C):" << endl;
	cin >> C(0) >> C(1) >> C(2);

	cout << "start point (A):\n" << A << endl;
	cout << "end point (C):\n" << C << endl;*/
	
	RowVector3d A(0.3, 0.15, 0.32);
	RowVector3d C(0.3, 0.15, 0.89);

	//Calculate the height of A and C
	d_xyz = inter.cal_dz(A, C);
	cout << "d_xyz:\n" << d_xyz << endl << endl;
	
	//Calculate point B
	B = inter.cal_pointB(A, C, d_xyz);
	cout << "B:\n" << B << endl << endl;

	//Calculate step size
	stepL = inter.cal_stepL(A, C);
	cout << "stepL:\n" << stepL << endl << endl;

	//Calculate center and radius
	center = inter.cal_center(A, B, C, &rad);
	cout << "center:\n" << center << endl;
	cout << "rad:\n" << rad << endl << endl;

	//Calculate the inter points
	points = inter.cal_Inter(A, B, C, center, rad, stepL);
	for (i = 0; i < points.cols(); i++)
	{
		cout << "ponit " << i + 1 << " : ";
		for (j = 0; j < 3; j++)
			cout << points(j, i) << " ";
		cout << endl;
	}

	cout << endl;

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "It costs: " << duration << "seconds\n";
	
	//避免运行结果框闪退
	system("pause");

	return 0;
}