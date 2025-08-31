// Angles header

#include <cmath>

class Angle {
	double r;
public:
	explicit Angle(double rad);
	double get() const;
	static Angle radians(double rad);
	static Angle degrees(double deg);
	double toRadians() const;
	double toDegrees() const;
};