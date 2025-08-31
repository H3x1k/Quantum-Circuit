// Angles implementation

#include "Angle.hpp"

#define PI   3.1415926535
#define DTR  0.0174532925
#define RTD 57.2957795130

Angle::Angle(double rad) : r(rad) {}

double Angle::get() const {
	return r;
}

Angle Angle::radians(double rad) {
	return Angle(rad);
}

Angle Angle::degrees(double deg) {
	return Angle(deg * DTR);
}

double Angle::toRadians() const {
	return r;
}

double Angle::toDegrees() const {
	return r * RTD;
}