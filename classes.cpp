#include <GraphicsMagick/Magick++.h>
#include <CImg.h>
#include "classes.h"

/** Vector Implementation. **/
Vector::Vector(float x, float y, float z) {
	_x = x;
	_y = y;
	_z = z;
}

Vector & Vector::operator+(Vector &v) {
	v._x += _x;
	v._y += _y;
	v._z += _z;
	return v;
}