#include <iostream>
#include <cmath>
using namespace std;
const double eps = 1e-10; //极小值

//取符号
int dcmp(double x){
    if (fabs(x) < eps) return 0;
    return x < 0 ? -1 : 1;
}

//点定义
struct Point
{
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};

//向量定义
typedef Point Vector;
Vector operator+(Vector A, Vector B){
    return Vector(A.x + B.x, A.y + B.y);
}
Vector operator-(Vector A, Vector B){
    return Vector(A.x - B.x, A.y - B.y);
}
Vector operator*(Vector A, double p){
    return Vector(A.x * p, A.y * p);
}
Vector operator/(Vector A, double p){
    return Vector(A.x / p, A.y / p);
}
bool operator<(Vector A, Vector B){
    return A.x < B.x || (A.x == B.x) && A.y < B.y;
}
bool operator==(Vector A, Vector B){
    return dcmp(A.x - B.x) == 0 && dcmp(A.y-B.y) == 0;
}

// !基本运算
double Dot(Vector A,Vector B){ //点乘
    return A.x * B.x + A.y * B.y;
}
double Length(Vector A){ //求长度
    return sqrt(Dot(A,A));
}
double Angle(Vector A, Vector B){ //求角度
    return acos(Dot(A,B))/Length(A)/Length(B);
}
double Cross(Vector A, Vector B){ //点乘
    return A.x * B.y - A.y * B.x;
}
double Area2(Point A, Point B, Point C){ //求三角面积
    return Cross(B-A,C-A);
}

// !点与直线
Point GetLineIntersection(Point P, Vector v, Point Q, Vector w){ 
    Vector u = P-Q;
    double t = Cross(w,u)/Cross(v,w);
    return P+v*t;
}//求直线交点
double DistanceToLine(Point P, Point A, Point B){
    Vector v1 = B - A,v2 = P - A;
    return fabs(Cross(v1, v2)) / Length(v1);
}//点到直线距离

// todo

// !凸包

//输入点数组p，输出为ch。函数返回凸包定点数
//输入不能有重复点
int ConvexHull(Point* p, int n, Point *ch){
    sort(p,p+n);
    int mm = 0;
    for (int i = 0;i < n;i ++){
        while(mm > 1 && Cross(ch[mm-1]-ch[mm-2],p[i]-ch[mm-2]) <= 0)
            mm--;
        ch[mm++] = p[i];
    }
    int k = mm;
    for (int i = n - 2;i >= 0;i --){
        while (mm > k && Cross(ch[mm-1] - ch[mm-2],p[i]-ch[mm-2]) <= 0)
            mm --;
        ch[mm++] = p[i];
    }
    if(n > 1) mm --;
    return mm;
}

int main()
{
}