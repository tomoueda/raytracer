#include <GraphicsMagick/Magick++.h>
#include <stdlib.h>
#include <vector>
#include <CImg.h>
#include <glog/logging.h>
#include "classes.h"

#define sqr(x) x*x
#define EPSILON .0001

/** Vector Implementation. &v is the actual object thus
    . is allowed, but actual v is a pointer in which case
    -> is needed to access elements. **/
Vector::Vector() {
	_x = 0;
	_y = 0;
	_z = 0;
}

Vector::Vector(float x, float y, float z) {
	_x = x;
	_y = y;
	_z = z;
}

/** Element-wise operations. **/

Vector & Vector::operator+(Vector &v) {
	float x = _x + v._x;
	float y = _y + v._y;
	float z = _z + v._z;
	return *new Vector(x, y, z);
}

void Vector::operator+=(Vector &v) {
	_x += v._x;
	_y += v._y;
	_z += v._z;
}

Vector & Vector::operator-(Vector &v) {
	float x = _x - v._x;
	float y = _y - v._y;
	float z = _z - v._z;
	return *new Vector(x, y, z);
}

void Vector::operator-=(Vector &v) {
	_x -= v._x;
	_y -= v._y;
	_z -= v._z;
}

Vector & Vector::operator*(float scalar) {
	float x = _x * scalar;
	float y = _y * scalar;
	float z = _z * scalar;
	return *new Vector(x, y, z);
}

void Vector::operator*=(float scalar) {
	_x *= scalar;
	_y *= scalar;
	_z *= scalar;
}


Vector & Vector::operator/(float scalar) {
	float x = _x / scalar;
	float y = _y / scalar;
	float z = _z / scalar;
	return *new Vector(x, y, z);
}

void Vector::operator/=(float scalar) {
	_x /= scalar;
	_y /= scalar;
	_z /= scalar;
}

bool Vector::operator==(const Vector &v) const {
	return fabs(_x - v._x) < EPSILON 
	&& fabs(_y - v._y) < EPSILON
	&& fabs(_z - v._z) < EPSILON;
}

bool Vector::operator!=(const Vector &v) const {
	return !(*this == v);
}

Normal & Vector::normalize() {
	Normal* n = new Normal(_x, _y, _z);
	return *n;
}

void Vector::setValue(float x, float y, float z) {
	_x = x;
	_y = y;
	_z = z;
}

void Vector::debug() {
	LOG(INFO) << "(" << _x << ", " << _y << ", " << _z << ", " << ")";
}

/** END VECTOR **/

/** 
	Begin Normal Implementation
**/

Normal::Normal() {
	_x = 0;
	_y = 0;
	_z = 0;
}
Normal::Normal(float x, float y, float z) {
	float magnitude = sqrt(sqr(x) + sqr(y) + sqr(z));
	if (magnitude > EPSILON) {
		_x = x / magnitude;
		_y = y / magnitude;
		_z = z / magnitude;
	} else {
		_x = 0;
		_y = 0;
		_z = 0;
	}
}

void Normal::normalize() {
	float magnitude = sqrt(sqr(_x) + sqr(_y) + sqr(_z));
	if (magnitude > EPSILON) {
		_x = _x / magnitude;
		_y = _y / magnitude;
		_z = _z / magnitude;
	} else {
		_x = 0;
		_y = 0;
		_z = 0;
	}
}

Normal& Normal::operator+(Normal& n) {
	float x = _x + n._x;
	float y = _y + n._y;
	float z = _z + n._z;
    this->normalize();
    return *new Normal(x, y, z);	
}

Normal& Normal::operator-(Normal& n) {
	float x = _x - n._x;
	float y = _y - n._y;
	float z = _z - n._z;
	this->normalize();
	return *new Normal(x, y, z);
}

void Normal::operator+=(Normal& n) {
	_x += n._x;
	_y += n._y;
	_z += n._z;
    this->normalize();
}

void Normal::operator-=(Normal& n) {
	_x -= n._x;
	_y -= n._y;
	_z -= n._z;
    this->normalize();
}

bool Normal::operator==(const Normal& n) const {
	return fabs(_x - n._x) < EPSILON 
	&& fabs(_y - n._y) < EPSILON
	&& fabs(_z - n._z) < EPSILON;
}

bool Normal::operator!=(const Normal& n) const {
	return !(*this == n);
}

void Normal::debug() {
	LOG(INFO) << "(" << _x << ", " << _y << ", " << _z << ", " << ")";
}

/** End Normal implementation **/

/** Start Point implementation. **/
Point::Point() {
	_x = 0;
	_y = 0;
	_z = 0;
}

Point::Point(float x, float y, float z) {
	_x = x;
	_y = y;
	_z = z;
}

Point& Point::operator+(Vector& v) {
	float x = _x + v._x;
	float y = _y + v._y;
	float z = _z + v._z;
	return *new Point(x, y, z);
}

Point& Point::operator-(Vector& v) {
	float x = _x - v._x;
	float y = _y - v._y;
	float z = _z - v._z;
	return *new Point(x, y, z);
}

void Point::operator+=(Vector& v) {
	_x += v._x;
	_y += v._y;
	_z += v._z;
}

void Point::operator-=(Vector& v) {
	_x -= v._x;
	_y -= v._y;
	_z -= v._z;
}

Vector& Point::operator-(Point& p) {
	float x = _x - p._x;
	float y = _y - p._y;
	float z = _z - p._z;
	return *new Vector(x, y, z);
}

bool Point::operator==(const Point& p) const {
	return fabs(_x - p._x) < EPSILON 
	&& fabs(_y - p._y) < EPSILON
	&& fabs(_z - p._z) < EPSILON;
}

bool Point::operator!=(const Point& p) const {
	return !(*this == p);
}

void Point::debug() {
	LOG(INFO) << "(" << _x << ", " << _y << ", " << _z << ", " << ")";
}

/** END IMPLMENTATION FOR POINT. **/


/**BEGIN IMPLEMENTATION FOR MATRIX + TRANSFORMATION. **/
/**END IMPLEMENTATION FOR MATRIX + TRANSFORMATION. **/


/** BEGIN IMPLEMENTATION FOR RAY. **/
Ray::Ray(Point &pos, Vector &dir, float t_min, float t_max) {
	_pos = pos;
	_dir = dir;
	_t_min = t_min;
	_t_max = t_max;
}

float Ray::get_t_min() {
	return _t_min;
}

float Ray::get_t_max() {
	return _t_max;
}

/** END IMPLEMENTATION FOR RAY. **/

/**BEGIN IMPLEMENTATION FOR MATRIX + TRANSFORMATION. **/
/**END IMPLEMENTATION FOR MATRIX + TRANSFORMATION. **/

/**BEGIN IMPLEMENTATION FOR COLOR. **/

Color::Color() {
	_r = 0;
	_b = 0;
	_g = 0;
}

Color::Color(float r, float g, float b) {
	_r = r;
	_b = b;
	_g = g;
}

/** Element-wise operations. **/

Color & Color::operator+(Color &c) {
	float r = _r + c._r;
	float g = _g + c._g;
	float b = _b + c._b;
	return *new Color(r, g, b);
}

Color & Color::operator-(Color &c) {
	float r = _r - c._r;
	float b = _b - c._b;
	float g = _g - c._g;
	return *new Color(r, g, b);
}

Color & Color::operator*(Color &c) {
	float r = _r * c._r;
	float b = _b * c._b;
	float g = _g * c._g;
	return *new Color(r, g, b);
}

Color & Color::operator/(Color &c) {
	float r = _r / c._r;
	float b = _b / c._b;
	float g = _g / c._g;
	return *new Color(r, g, b);
}

Color & Color::operator*(float scalar) {
	float r = _r * scalar;
	float b = _b * scalar;
	float g = _g * scalar;
	return *new Color(r, g, b);
}

Color & Color::operator/(float scalar) {
	float r = _r / scalar;
	float b = _b / scalar;
	float g = _g / scalar;
	return *new Color(r, g, b);
}

void Color::operator+=(Color &c) {
	_r += c._r;
	_g += c._g;
	_b += c._b;
}

void Color::operator-=(Color &c) {
	_r -= c._r;
	_g -= c._g;
	_b -= c._b;
}

void Color::operator*=(Color &c) {
	_r *= c._r;
	_g *= c._g;
	_b *= c._b;
}

void Color::operator/=(Color &c) {
	_r /= c._r;
	_g /= c._g;
	_b /= c._b;
}

void Color::operator*=(float scalar) {
	_r *= scalar;
	_g *= scalar;
	_b *= scalar;
}

void Color::operator/=(float scalar) {
	_r /= scalar;
	_g /= scalar;
	_b /= scalar;
}

bool Color::operator==(const Color &c) const {
	return fabs(_r - c._r) < EPSILON 
	&& fabs(_g - c._g) < EPSILON
	&& fabs(_b - c._b) < EPSILON;
}

bool Color::operator!=(const Color &c) const {
	return !(*this == c);
}

void Color::setValue(float r, float g, float b) {
	_r = r;
	_g = g;
	_b = b;
}

void Color::debug() {
	LOG(INFO) << "(" << _r << ", " << _g << ", " << _b << ", " << ")";
}

/** Takes in a vector V and takes the values of v and
	changes it to r, g, b values of this color. **/
void Color::convert(Vector &v) {
	_r = v._x;
	_g = v._y;
	_b = v._z;
}


/**END IMPLEMENTATION FOR COLOR. **/

/**BEGIN IMPLEMENTATION FOR BRDF. **/
BRDF::BRDF(Color& kd, Color& ks, Color& ka, Color& kr) {
	_kd = kd;
	_ks = ks;
	_ka = ka;
	_kr = kr;
}

Color& BRDF::get_kd() {
	return _kd;
}

Color& BRDF::get_ks() {
	return _ks;
}

Color& BRDF::get_ka() {
	return _ka;
}

Color& BRDF::get_kr() {
	return _kr;
}

/**END IMPLEMENTATION FOR BRDF. **/

/**BEGIN IMPLEMENTATION FOR SAMPLE. **/
Sample::Sample(float x, float y) {
	_x = x;
	_y = y;
}

float Sample::get_x() {
	return _x;
}

float Sample::get_y() {
	return _y;
}
/**END IMPLEMENTATION FOR SAMPLE. **/

/**BEGIN IMPLEMENTATION FOR LOCAL GEO. **/
LocalGeo::LocalGeo(Point &pos, Point &normal) {
	_pos = pos;
	_normal = normal;
}

Point& LocalGeo::get_pos() {
	return _pos;
}

Normal& LocalGeo::get_norm() {
	return _normal;
}
/**END IMPLEMENTATION FOR LOCAL GEO. **/

/** BEGIN IMPLEMENTATION FOR INTERSECTION. **/
Intersection::Intersection(LocalGeo& localGeo, Primitive* primitives) {
	_localGeo = localGeo;
	_primitive = primitive;
}

LocalGeo& Intersection::get_localGeo() {
	return _localGeo
}

Primitive* Intersection::get_primitive() {
	return _primitves;
}
/** END IMPLEMENTATION FOR INTERSECTION. **/

/**BEGIN IMPLEMENTATION FOR GEOMETRIC PRIMITIVE. **/
/**END IMPLEMENTATION FOR GEOMETRIC PRIMITIVE. **/

/**BEGIN IMPLEMENTATION FOR AGGREGATE PRIMITIVE. **/
/**END IMPLEMENTATION FOR AGGREGATE PRIMITIVE. **/

/**BEGIN IMPLEMENTATION FOR MATERIAL. **/
/**END IMPLEMENTATION FOR MATERIAL. **/

/**BEGIN IMPLEMENTATION FOR SAMPLER. **/
/**END IMPLEMENTATION FOR SAMPLER. **/

/**BEGIN IMPLEMENTATION FOR CAMERA. **/
/**END IMPLEMENTATION FOR CAMERA. **/

/**BEGIN IMPLEMENTATION FOR RAY TRACER. **/
/**END IMPLEMENTATION FOR RAY TRACER. **/

/**BEGIN IMPLEMENTATION FOR LIGHTS. **/
/**END IMPLEMENTATION FOR LIGHTS. **/

/**BEGIN IMPLEMENTATION FOR FILM. **/
/**END IMPLEMENTATION FOR FILM. **/

/**BEGIN IMPLEMENTATION FOR SCENE.  **/
/**END IMPLEMENTATION FOR SCENE. **/

