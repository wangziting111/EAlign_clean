#include "XYZ.h"

//-------------- constructor --------------//
XYZ::XYZ(void)
{
	X = 0;
	Y = 0;
	Z = 0;
}
XYZ::~XYZ(void)
{
}
XYZ::XYZ(double w)
{
	X = w;
	Y = w;
	Z = w;
}
XYZ::XYZ(double x, double y, double z)
{
	X = x;
	Y = y;
	Z = z;
}
XYZ::XYZ(const XYZ & xyz)
{
	this->X = xyz.X;
	this->Y = xyz.Y;
	this->Z = xyz.Z;
}
XYZ & XYZ::operator =(const XYZ & xyz)
{
	if(this==&xyz)return *this;
	this->X = xyz.X;
	this->Y = xyz.Y;
	this->Z = xyz.Z;
	return *this;
}

//---------------operator_reload ----------------------// 
XYZ XYZ::operator +(const XYZ & xyz) const
{
	XYZ temp;
	temp.X = X + xyz.X;
	temp.Y = Y + xyz.Y;
	temp.Z = Z + xyz.Z;
	return (temp);
}
XYZ XYZ::operator -(const XYZ & xyz) const
{
	XYZ temp;
	temp.X = X - xyz.X;
	temp.Y = Y - xyz.Y;
	temp.Z = Z - xyz.Z;
	return (temp);
}
XYZ XYZ::operator *(double i) const
{
	XYZ temp;
	temp.X = X * i;
	temp.Y = Y * i;
	temp.Z = Z * i;
	return (temp);
}
XYZ XYZ::operator /(double i) const
{
	XYZ temp;
	temp.X = X / i;
	temp.Y = Y / i;
	temp.Z = Z / i;
	return (temp);
}
void XYZ::operator +=(const XYZ & xyz)
{
	this->X += xyz.X;
	this->Y += xyz.Y;
	this->Z += xyz.Z;
}
void XYZ::operator -=(const XYZ & xyz)
{
	this->X -= xyz.X;
	this->Y -= xyz.Y;
	this->Z -= xyz.Z;
}
void XYZ::operator *=(double i)
{
	this->X *= i;
	this->Y *= i;
	this->Z *= i;
}
void XYZ::operator /=(double i)
{
	this->X /= i;
	this->Y /= i;
	this->Z /= i;
}
void XYZ::operator =(double i)
{
	this->X = i;
	this->Y = i;
	this->Z = i;
}

//---------------function ----------------------// 
void XYZ::XYZ2()
{
	X = X * X;
	Y = Y * Y;
	Z = Z * Z;
}
double XYZ::distance(const XYZ & xyz) const
{
	XYZ xyz_ = xyz;
	xyz_ -= *this;
	xyz_.XYZ2();
	double sum = xyz_.X + xyz_.Y + xyz_.Z;
	return (sum > 0.0 ? sqrt(sum) : 0.0);
}
double XYZ::distance_square(const XYZ & xyz)  const 
{
	XYZ xyz_ = xyz;
	xyz_ -= *this;
	xyz_.XYZ2();
	double sum = xyz_.X + xyz_.Y + xyz_.Z;
	return (sum);
}
void XYZ::xyz2double(double t[3]) const 
{
	t[0] = this->X;
	t[1] = this->Y;
	t[2] = this->Z;
} 
void XYZ::double2xyz(double t[3])
{
	this->X = t[0];
	this->Y = t[1];
	this->Z = t[2];
}

//---- angle related ----//

