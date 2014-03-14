// #include <GraphicsMagick/Magick++.h>
#include <stdlib.h>
#include <vector>
#include <CImg.h>
#include <glog/logging.h>
#undef Success
#include <Eigen/Dense>
#include "classes.h"

#define EPSILON .0001
using namespace Eigen;
using namespace cimg_library;

/** STATIC FUNCTIONS **/
inline float sqr(float x) {
    return x * x;
}

inline float max(float a, float b) {
    if (a > b) return a;
    return b;
}

inline float max(float a, float b, float c) {
    if (a > b and a > c) {
        return a;
    }
    if (b > c) {
        return b;
    }
    return c;
}

inline float min(float a, float b) {
    if (a < b) return a;
    return b;
}

inline float min(float a, float b, float c) {
    if (a < b and a < c) {
        return a;
    }
    if (b < c) {
        return b;
    }
    return c;
}

CImg<float> scaleForJPG(CImg<float> img) {
    return ((img - img.max()) / (img.max() - img.min())) * 255;
}


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

//Used for shading
Vector & Vector::operator+(Normal &v) {
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

Vector& Vector::cross_product(Vector& v) {
    float x = _y * v._z - _z * v._y;
    float y = _z * v._x - _x * v._z;
    float z = _x * v._y - _y * v._x;
    return *new Vector(x, y, z);
}

float Vector::get_x() {
    return _x;
}

float Vector::get_y() {
    return _y;
}

float Vector::get_z() {
    return _z;
}

Normal & Vector::get_normal() {
    Normal* n = new Normal(_x, _y, _z);
    return *n;
}

void Vector::normalize() {
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

void Vector::setValue(float x, float y, float z) {
    _x = x;
    _y = y;
    _z = z;
}

void Vector::setValue(Normal& n) {
    _x = n._x;
    _y = n._y;
    _z = n._z;
}

float Vector::dot_product(Vector& v) {
    return _x * v._x + _y * v._y + _z * v._z;
}

float Vector::dot_product(Normal& v) {
    return _x * v._x + _y * v._y + _z * v._z;
}

void Vector::debug() {
    LOG(INFO) << "(" << _x << ", " << _y << ", " << _z << ")";
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

/** Again convenient for shading. **/
Normal& Normal::operator*(float s) {
    _x *= s;
    _y *= s;
    _z *= s;
    this->normalize();
    float x = _x;
    float y = _y;
    float z = _z;
    return * new Normal(x, y, z);
}

bool Normal::operator==(const Normal& n) const {
    return fabs(_x - n._x) < EPSILON 
    && fabs(_y - n._y) < EPSILON
    && fabs(_z - n._z) < EPSILON;
}

bool Normal::operator!=(const Normal& n) const {
    return !(*this == n);
}
 
float Normal::get_x() {
    return _x;
}

float Normal::get_y() {
    return _y;
}

float Normal::get_z() {
    return _z;
}

float Normal::dot_product(Normal &n) {
    return _x * n._x + _y * n._y + _z * n._z;
}

float Normal::dot_product(Vector &v) {
    return _x * v._x + _y * v._y + _z * v._z;
}

void Normal::debug() {
    LOG(INFO) << "(" << _x << ", " << _y << ", " << _z << ")";
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

/** This is technically nonsense, but is useful for finding p
    in the camera. **/
Point& Point::operator*(float f) {
    float x = _x * f;
    float y = _y * f;
    float z = _z * f;
    return *new Point(x, y, z);
}

/** This is technically nonsense, but is useful for finding p
    in the camera. **/
Point& Point::operator+(Point& p) {
    float x = _x + p._x;
    float y = _y + p._y;
    float z = _z + p._z;
    return *new Point(x, y, z);
}

bool Point::operator==(const Point& p) const {
    return fabs(_x - p._x) < EPSILON 
    && fabs(_y - p._y) < EPSILON
    && fabs(_z - p._z) < EPSILON;
}

bool Point::operator!=(const Point& p) const {
    return !(*this == p);
}

float Point::get_x() {
    return _x;
}

float Point::get_y() {
    return _y;
}

float Point::get_z() {
    return _z;
}

float Point::dot_product(Point& p) {
    return _x * p._x + _y * p._y + _z * p._z;
}

float Point::dot_product(Normal& n) {
    return _x * n._x + _y * n._y + _z * n._z;
}

float Point::distance(Point& p) {
    float ans = sqrt(sqr(_x - p._x) + sqr(_y - p._y) + sqr(_z - p._z));
    return ans;
}

void Point::debug() {
    LOG(INFO) << "(" << _x << ", " << _y << ", " << _z << ")";
}

/** END IMPLMENTATION FOR POINT. **/


/**BEGIN IMPLEMENTATION FOR MATRIX + TRANSFORMATION. **/

/** This constructor will create an identity matrix.
    If applied this matrix will give a identity transformation.
    REMEMBER THAT THE MATRIX TENDS TO EXECUTE TRANSFORMATION IN THE
    ORDER THEY ARE PUT IN. **/
TransMatrix::TransMatrix() {
    _mat << 1, 0, 0, 0,
           0, 1, 0, 0,
           0, 0, 1, 0,
           0, 0, 0, 1;
}

/** A constructor that takes in a MAT which is a Matrix4f and copies it
    to our member. **/
TransMatrix::TransMatrix(Matrix4f mat) {
    _mat = mat;
}

/** Adds a rotation to the current matrix. This matrix will rotate points
    about ABOUT by radian RADIAN. Remember that quaternion rotation matrix
    are NOT commutative, so take care when deciding the order of rotation. **/
void TransMatrix::add_rotation(Point& about, float radian) {
    Matrix4f new_rotate;
    Matrix4f rx;
    Matrix4f identity;
    float s = sin(radian);
    float c = cos(radian);
    float x = about.get_x();
    float y = about.get_y();
    float z = about.get_z();
    identity << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
    rx << 0, -1*z, y, 0,
          z, 0, -1*x, 0,
          -1*y, x, 0, 0,
          0, 0, 0, 1;
    new_rotate = identity + (rx)*s + (rx)*(rx)*(1-c);
    // new_rotate << c+x*x*(1-c), x*y*(1-c)-z*s, x*z*(1-c)+y*s, 0,
    //               y*x*(1-c)+z*s, c+y*y*(1-c), y*z*(1-c)-x*s, 0,
    //               z*x*(1-c)-y*s, z*y*(1-c)+x*s, c+z*z*(1-c), 0,
    //               0, 0, 0, 1;
    _mat = _mat * new_rotate;
}

/** A translation matrix that will move the point by, X, Y, Z. **/
void TransMatrix::add_translation(float x, float y, float z) {
    Matrix4f new_translate;
    new_translate << 1, 0, 0, x,
                     0, 1, 0, y,
                     0, 0, 1, z, 
                     0, 0, 0, 1;
    _mat = _mat * new_translate;
}

/** Adds on a transformation that will scale points by X, Y, Z. **/
void TransMatrix::add_scaling(float x, float y, float z) {
    Matrix4f new_scaling;
    new_scaling << x, 0, 0, 0,
                   0, y, 0, 0,
                   0, 0, z, 0,
                   0, 0, 0, 1;
    _mat = _mat * new_scaling;
}

TransMatrix& TransMatrix::invt() {
    return *new TransMatrix(_mat.inverse().transpose());
}

TransMatrix& TransMatrix::inv() {
    return *new TransMatrix(_mat.inverse());
}

Vector4f TransMatrix::operator*(Vector4f v) {
    return _mat * v;
}

Transformation::Transformation() {
    //Empty constructor.
}

/** This constructor takes in M and stores the original to transform
    points and also keeps a copy of the inverse transpose. **/
Transformation::Transformation(TransMatrix& m) {
    _m = m;
    _minvt = m.invt();
}

/** Overloaded Constructors. **/
Point& Transformation::operator*(Point& p) {
    Vector4f point(p.get_x(), p.get_y(), p.get_z(), 1);
    point = _m * point;
    return *new Point(point(0), point(1), point(2));
}

Vector& Transformation::operator*(Vector& v) {
    Vector4f vec(v.get_x(), v.get_y(), v.get_z(), 0);
    vec = _m * vec;
    return *new Vector(vec(0), vec(1), vec(2));
}

Normal& Transformation::operator*(Normal& n) {
    Vector4f norm(n.get_x(), n.get_y(), n.get_z(), 0);
    norm = _minvt * norm;
    return *new Normal(norm(0), norm(1), norm(2));
}

Ray& Transformation::operator*(Ray& r) {
    return *new Ray(*this * r.get_pos(), *this * r.get_dir(), 
        r.get_t_min(), r.get_t_max());
}

LocalGeo& Transformation::operator*(LocalGeo& l) {
    return *new LocalGeo(*this * l.get_pos(), *this * l.get_norm());
}

Transformation& Transformation::inv() {
    return *new Transformation(_m.inv());
}


/**END IMPLEMENTATION FOR MATRIX + TRANSFORMATION. **/


/** BEGIN IMPLEMENTATION FOR RAY. **/
Ray::Ray() {
    //Default Constructor
}

Ray::Ray(Point& pos, Vector& dir, float t_min, float t_max) {
    _pos = pos;
    _dir = dir;
    _dir.normalize();
    _t_min = t_min;
    _t_max = t_max;
}

float Ray::get_t_min() {
    return _t_min;
}

float Ray::get_t_max() {
    return _t_max;
}

Point& Ray::get_pos() {
    return _pos;
}

Vector& Ray::get_dir() {
    return _dir;
}

Point& Ray::get_pos_with_t(float t) {
    return _pos + _dir * t;
}

void Ray::update_tmax(float t) {
    _t_max = t;
}

void Ray::update_origin(Point& o) {
    _pos = o;
}

void Ray::update_dir(Vector& dir) {
    _dir = dir;
    _dir.normalize();
}

void Ray::update_tmin(float t) {
    _t_min = t;
}

/** END IMPLEMENTATION FOR RAY. **/

/**BEGIN IMPLEMENTATION FOR MATRIX + TRANSFORMATION. **/
/**END IMPLEMENTATION FOR MATRIX + TRANSFORMATION. **/

/**BEGIN IMPLEMENTATION FOR COLOR. **/

Color::Color() {
    _r = 0;
    _g = 0;
    _b = 0;
}

Color::Color(float r, float g, float b) {
    _r = r;
    _g = g;
    _b = b;
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

void Color::setValue(Color& c) {
    _r = c._r;
    _g = c._g;
    _b = c._b;
}

void Color::black() {
    _r = 0;
    _g = 0;
    _b = 0;
}

void Color::white() {
    _r = 254;
    _g = 254;
    _b = 254;
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

float Color::get_r() {
    return _r;
}

float Color::get_g() {
    return _g;
}

float Color::get_b() {
    return _b;
}

/**END IMPLEMENTATION FOR COLOR. **/

/** Begin implementation for vertex. **/
Vertex::Vertex() {
}

Vertex::Vertex(Point& pos, Normal& norm) {
    _pos = pos;
    _norm = norm;
}

Point& Vertex::get_point() {
    return _pos;
}

Normal& Vertex::get_norm() {
    return _norm;
}

void Vertex::set_pos(Point& pos) {
    _pos = pos;
}

void Vertex::set_norm(Normal& norm) {
    _norm = norm;
}
/** End implementation for vertex. **/

/**BEGIN IMPLEMENTATION FOR BRDF. **/
BRDF::BRDF() {
    //empty constructor.
}

BRDF::BRDF(Color& kd, Color& ks, Color& ka, Color& kr, float sp) {
    _kd = kd;
    _ks = ks;
    _ka = ka;
    _kr = kr;
    _sp = sp;
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

float BRDF::get_sp() {
    return _sp;
}

void BRDF::setValue(BRDF& brdf) {
    _kd = brdf._kd;
    _ks = brdf._ks;
    _ka = brdf._ka;
    _kr = brdf._kr;
    _sp = brdf._sp;
} 

bool BRDF::isReflect() {
    return _kr.get_r() > 0 && _kr.get_g() > 0 && _kr.get_b() > 0;
}

/**END IMPLEMENTATION FOR BRDF. **/

/**BEGIN IMPLEMENTATION FOR SAMPLE. **/
Sample::Sample() {
    //Default Constructor
}

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

void Sample::update_x(float x) {
    _x = x;
}

void Sample::update_y(float y) {
    _y = y;
}
/**END IMPLEMENTATION FOR SAMPLE. **/

/**BEGIN IMPLEMENTATION FOR LOCAL GEO. **/
LocalGeo::LocalGeo() {
    /** Empty Constructor. **/
}

LocalGeo::LocalGeo(Point &pos, Normal &normal) {
    _pos = pos;
    _normal = normal;
}

Point& LocalGeo::get_pos() {
    return _pos;
}

Normal& LocalGeo::get_norm() {
    return _normal;
}

void LocalGeo::set_pos(Point& p) {
    _pos = p;
}

void LocalGeo::set_norm(Normal& n) {
    _normal = n;
}
/**END IMPLEMENTATION FOR LOCAL GEO. **/

/** BEGIN SHAPE IMPLEMENTATION. **/

Sphere::Sphere(Point& center, float radius) {
    _center = center;
    _radius = radius;
}

/** The Intersect function for the Sphere. It takes in a RAY,
    which contain tmin and tmax, and the origin and direction of the ray.
    We also pass in THIT, a pointer to a float, which we will update
    the value if the ray hits this sphere. LOCAL is a pointer to a
    LocalGeo which we will update with Shading information if the ray
    intersects the Sphere. We are simply plugging in the function of the
    ray to the explicit formula of the sphere and using the quadratic formula
    to find intersections. The function itself will return TRUE if 
    the ray and this sphere interesects and false if it doesn't. **/
bool Sphere::intersect(Ray& ray, float *thit, LocalGeo* local) {
    float mint = ray.get_t_min();
    float maxt = ray.get_t_max();
    Point pos = ray.get_pos();
    Vector dir = ray.get_dir();
    Vector sphere_dir = pos - _center;
    float a = dir.dot_product(dir);
    float b = 2 * dir.dot_product(sphere_dir);
    float c = sphere_dir.dot_product(sphere_dir) - sqr(_radius);
    float determinant = sqr(b) - 4 * a * c;
    if (determinant < 0) {
        return false;
    }
    float tnew1 = (-b + sqrt(determinant)) / (2*a);
    float tnew2 = (-b - sqrt(determinant)) / (2*a);
    if (tnew1 <= tnew2) {
        if (tnew1 < mint || tnew2 > maxt) {
            return false;
        }
        *thit = tnew1;
    } else {
        if (tnew2 < mint || tnew1 > maxt) {
            return false;
        }
        *thit = tnew2;
    }
    Point intersect = ray.get_pos_with_t(*thit);
    Vector intersect_norm_dir = intersect - _center;
    Normal intersect_norm = intersect_norm_dir.get_normal();
    local->set_pos(intersect);
    local->set_norm(intersect_norm);
    return true;
}

/** This function is used mainly for shadow rays. The function takes in a 
    RAY just like the function above, but we don't need to update
    any of the variables. **/
bool Sphere::intersectP(Ray& ray) {
    float mint = ray.get_t_min();
    float maxt = ray.get_t_max();
    Point pos = ray.get_pos();
    Vector dir = ray.get_dir();
    Vector sphere_dir = pos - _center;
    float a = dir.dot_product(dir);
    float b = 2 * dir.dot_product(sphere_dir);
    float c = sphere_dir.dot_product(sphere_dir) - sqr(_radius);
    float determinant = sqr(b) - 4 * a * c;
    if (determinant < 0 || a == 0) {
        return false;
    }
    float tnew1 = (-b + sqrt(determinant)) / (2*a);
    float tnew2 = (-b - sqrt(determinant)) / (2*a);
    if (tnew1 <= tnew2) {
        if (tnew1 < mint || tnew1 > maxt) {
            return false;
        }
    } else {
        if (tnew2 < mint || tnew2 > maxt) {
            return false;
        }
    }
    return true;
}

BoundingBox& Sphere::createBoundingBox() {
    float xmax = _center.get_x() + _radius;
    float ymax = _center.get_y() + _radius;
    float zmax = _center.get_z() + _radius;
    float xmin = _center.get_x() - _radius;
    float ymin = _center.get_y() - _radius;
    float zmin = _center.get_z() - _radius;
    return *new BoundingBox(xmax, ymax, zmax, xmin, ymin, zmin);
}


Point& Sphere::getCenter() {
    return _center;
}

float Sphere::getRadius() {
    return _radius;
}

/** The constructor for triangle takes in V1, V2, V3, which are
    Points but we will represent them as vertices. **/
Triangle::Triangle(Point& v1, Point& v2, Point& v3) {
    _v1 = v1;
    _v2 = v2;
    _v3 = v3;
}

Point& Triangle::getv1() {
    return _v1;
}

Point& Triangle::getv2() {
    return _v2;
}

Point& Triangle::getv3() {
    return _v3;
}


/** This function takes in RAY a ray and returns true if the ray
    intersects this triangle and false if the ray misses.
    If the ray does hit this triangle we
    update THIT, with the intersection point, and we also update,
    LOCAL with the point of intersection and norm of the intersectin. **/
bool Triangle::intersect(Ray& ray, float* thit, LocalGeo* local) {
    float mint = ray.get_t_min();
    float maxt = ray.get_t_max();
    Vector g = ray.get_dir() * -1;
    Vector j = _v2 - _v1;
    Vector k = _v3 - _v1;
    Vector l = ray.get_pos() - _v1;
    float gx = g.get_x();
    float gy = g.get_y();
    float gz = g.get_z();
    float jx = j.get_x();
    float jy = j.get_y();
    float jz = j.get_z();
    float kx = k.get_x();
    float ky = k.get_y();
    float kz = k.get_z();
    float lx = l.get_x();
    float ly = l.get_y();
    float lz = l.get_z();
    Matrix3f A;
    A << jx, kx, gx,
         jy, ky, gy,
         jz, kz, gz;
    Vector3f b(lx, ly, lz);
    Vector3f x = A.inverse() * b;
    float gamma = x(0);
    float beta = x(1);
    float t = x(2);
    if (gamma + beta > 1 || gamma < 0 || beta < 0 || 
        t < mint || t > maxt) {
        return false;
    }
    *thit = t;
    Point pos = ray.get_pos_with_t(t);
    Vector norm_dir = k.cross_product(j);
    Normal norm = norm_dir.get_normal();
    local->set_norm(norm);
    local->set_pos(pos);
    return true;
}

/** This function takes in a RAY and determines if it intersects with
    this triangle. If it does intersect we return true. Unlike the previous
    function we don't update any parameters, because this is primarily for
    the shadow rays. **/
bool Triangle::intersectP(Ray& ray) {
    Vector nonmoving(0, 0, 0);
    Vector dir = ray.get_dir();
    if (dir == nonmoving) {
        return false;
    }
    float mint = ray.get_t_min();
    float maxt = ray.get_t_max();
    Vector g = dir * -1;
    Vector j = _v2 - _v1;
    Vector k = _v3 - _v1;
    Vector l = ray.get_pos() - _v1;
    float gx = g.get_x();
    float gy = g.get_y();
    float gz = g.get_z();
    float jx = j.get_x();
    float jy = j.get_y();
    float jz = j.get_z();
    float kx = k.get_x();
    float ky = k.get_y();
    float kz = k.get_z();
    float lx = l.get_x();
    float ly = l.get_y();
    float lz = l.get_z();
    Matrix3f A;
    A << jx, kx, gx,
         jy, ky, gy,
         jz, kz, gz;
    Vector3f b(lx, ly, lz);
    Vector3f x = A.inverse() * b;
    float gamma = x(0);
    float beta = x(1);
    float t = x(2);
    if (gamma + beta > 1 || gamma < 0 || beta < 0 || 
        t < mint || t > maxt) {
        return false;
    }
    return true;
}

BoundingBox& Triangle::createBoundingBox() {
    float xmax = max(_v1.get_x(), _v2.get_x(), _v3.get_x());
    float ymax = max(_v1.get_y(), _v2.get_y(), _v3.get_y());
    float zmax = max(_v1.get_z(), _v2.get_z(), _v3.get_z());
    float xmin = min(_v1.get_x(), _v2.get_x(), _v3.get_x());
    float ymin = min(_v1.get_y(), _v2.get_y(), _v3.get_y());
    float zmin = min(_v1.get_z(), _v2.get_z(), _v3.get_z());
    return *new BoundingBox(xmax, ymax, zmax, xmin, ymin, zmin);
}


VertexTriangle::VertexTriangle(Vertex& v1, Vertex& v2, Vertex& v3) {
    _v1 = v1;
    _v2 = v2;
    _v3 = v3;
}

bool VertexTriangle::intersect(Ray& ray, float* thit, LocalGeo* local) {
    float mint = ray.get_t_min();
    float maxt = ray.get_t_max();
    Vector g = ray.get_dir() * -1;
    Vector j = _v2.get_point() - _v1.get_point();
    Vector k = _v3.get_point() - _v1.get_point();
    Vector l = ray.get_pos() - _v1.get_point();
    float gx = g.get_x();
    float gy = g.get_y();
    float gz = g.get_z();
    float jx = j.get_x();
    float jy = j.get_y();
    float jz = j.get_z();
    float kx = k.get_x();
    float ky = k.get_y();
    float kz = k.get_z();
    float lx = l.get_x();
    float ly = l.get_y();
    float lz = l.get_z();
    Matrix3f A;
    A << jx, kx, gx,
         jy, ky, gy,
         jz, kz, gz;
    Vector3f b(lx, ly, lz);
    Vector3f x = A.inverse() * b;
    float gamma = x(0);
    float beta = x(1);
    float t = x(2);
    if (gamma + beta > 1 || gamma < 0 || beta < 0 || 
        t < mint || t > maxt) {
        return false;
    }
    *thit = t;
    Point pos = ray.get_pos_with_t(t);
    float alpha = 1 - (gamma + beta);
    Normal norm = (_v1.get_norm() * alpha) + (_v2.get_norm() * beta) + (_v3.get_norm() * gamma);
    local->set_norm(norm);
    local->set_pos(pos);
    return true;
}

bool VertexTriangle::intersectP(Ray& ray) {
    float mint = ray.get_t_min();
    float maxt = ray.get_t_max();
    Vector g = ray.get_dir() * -1;
    Vector j = _v2.get_point() - _v1.get_point();
    Vector k = _v3.get_point() - _v1.get_point();
    Vector l = ray.get_pos() - _v1.get_point();
    float gx = g.get_x();
    float gy = g.get_y();
    float gz = g.get_z();
    float jx = j.get_x();
    float jy = j.get_y();
    float jz = j.get_z();
    float kx = k.get_x();
    float ky = k.get_y();
    float kz = k.get_z();
    float lx = l.get_x();
    float ly = l.get_y();
    float lz = l.get_z();
    Matrix3f A;
    A << jx, kx, gx,
         jy, ky, gy,
         jz, kz, gz;
    Vector3f b(lx, ly, lz);
    Vector3f x = A.inverse() * b;
    float gamma = x(0);
    float beta = x(1);
    float t = x(2);
    if (gamma + beta > 1 || gamma < 0 || beta < 0 || 
        t < mint || t > maxt) {
        return false;
    }
    return true;
}

Vertex& VertexTriangle::getv1() {
    return _v1;
}

Vertex& VertexTriangle::getv2() {
    return _v2;
}

Vertex& VertexTriangle::getv3() {
    return _v3;
}

BoundingBox& VertexTriangle::createBoundingBox() {
    float xmax = max(_v1.get_point().get_x(),
       _v2.get_point().get_x(), _v3.get_point().get_x());
    float ymax = max(_v1.get_point().get_y(),
       _v2.get_point().get_y(), _v3.get_point().get_y());
    float zmax = max(_v1.get_point().get_z(),
       _v2.get_point().get_z(), _v3.get_point().get_z());
    float xmin = min(_v1.get_point().get_x(),
       _v2.get_point().get_x(), _v3.get_point().get_x());
    float ymin = min(_v1.get_point().get_y(),
       _v2.get_point().get_y(), _v3.get_point().get_y());
    float zmin = min(_v1.get_point().get_z(),
       _v2.get_point().get_z(), _v3.get_point().get_z());
    return *new BoundingBox(xmax, ymax, zmax, xmin, ymin, zmin);
}


BoundingBox::BoundingBox() {
    //Default Constructor.
}

BoundingBox::BoundingBox(BoundingBox& bb1, BoundingBox& bb2) {
    _xmax = max(bb1._xmax, bb2._xmax);
    _ymax = max(bb1._ymax, bb2._ymax);
    _zmax = max(bb1._zmax, bb2._zmax);
    _xmin = min(bb1._xmin, bb2._xmin);
    _ymin = min(bb1._ymin, bb2._ymin);
    _zmin = min(bb1._zmin, bb2._zmin);
}

BoundingBox::BoundingBox(float xmax, float ymax, float zmax,
    float xmin, float ymin, float zmin) {
    _xmax = xmax;
    _ymax = ymax;
    _zmax = zmax;
    _xmin = xmin;
    _ymin = ymin;
    _zmin = zmin;
}

bool BoundingBox::intersect(Ray& ray, float* thit, LocalGeo* local) {
    exit(1);
    //This should never be called. 
}

bool BoundingBox::intersectP(Ray& ray) {
    Point rayPos = ray.get_pos();
    Vector rayDir = ray.get_dir();
    float t_xmin, t_xmax, t_ymin, t_ymax, t_zmin, t_zmax;
    float alphax = 1 / rayDir.get_x();
    if (alphax >= 0) {
        t_xmin = alphax * (_xmin - rayPos.get_x());
        t_xmax = alphax * (_xmax - rayPos.get_x());
    } else {
        t_xmin = alphax * (_xmax - rayPos.get_x());
        t_xmax = alphax * (_xmin - rayPos.get_x());
    }
    float alphay = 1 / rayDir.get_y();
    if (alphay >= 0) {
        t_ymin = alphay * (_ymin - rayPos.get_y());
        t_ymax = alphay * (_ymax - rayPos.get_y());
    } else {
        t_ymin = alphay * (_ymax - rayPos.get_y());
        t_ymax = alphay * (_ymin - rayPos.get_y());
    }

    float alphaz = 1 / rayDir.get_z();
    if (alphaz >= 0) {
        t_zmin = alphaz * (_zmin - rayPos.get_z());
        t_zmax = alphaz * (_zmax - rayPos.get_z());
    } else {
        t_zmin = alphaz * (_zmax - rayPos.get_z());
        t_zmax = alphaz * (_zmin - rayPos.get_z());
    }
    if (t_xmin > t_ymax || t_ymin > t_zmax || t_zmin > t_xmax) {
        return false;
    }
    return true;
}

BoundingBox& BoundingBox::createBoundingBox() {
    return *this;
}

/** Return the center of the boundingbox along AXIS.
    Axis = 0 corresponds to the x-axis
    Axis = 1 corresponds to the y-axis
    Axis = 2 corresponds to the z-axis.
    Returns null if anything else if input. **/
float BoundingBox::getCenter(int axis) {
    if (axis == 0) {
        return (_xmax + _xmin) / 2;
    } else if (axis == 1) {
        return (_ymax + _ymin) / 2;
    } else if (axis == 2) {
        return (_zmax + _zmin) / 2;
    }
    return (float) NULL;
}

/** END SHAPE IMPLEMENTATION. **/


/** BEGIN IMPLEMENTATION FOR INTERSECTION. **/
Intersection::Intersection() {
    //Default Constructor
}

Intersection::Intersection(LocalGeo& localGeo, Primitive* primitive) {
    _localGeo = localGeo;
    _primitive = primitive;
}

LocalGeo& Intersection::get_localGeo() {
    return _localGeo;
}

Primitive* Intersection::get_primitive() {
    return _primitive;
}

void Intersection::set_primitive(Primitive* primitive) {
    _primitive = primitive;
}

void Intersection::set_local(LocalGeo& localGeo) {
    _localGeo = localGeo;
}
/** END IMPLEMENTATION FOR INTERSECTION. **/

/**BEGIN IMPLEMENTATION FOR GEOMETRIC PRIMITIVE. **/
GeometricPrimitive::GeometricPrimitive() {
    //Empty constructor
}

/** This constructor takes in the object's transformation, TRANS, and the
    shape of the primitive, and the material of the object. **/
GeometricPrimitive::GeometricPrimitive(Transformation& trans, Shape* shape, Material* mat) {
    _objToWorld = trans;
    _worldToObj = trans.inv();
    _shape = shape;
    _mat = mat;
}

bool GeometricPrimitive::intersect(Ray& ray, float *thit, Intersection* in) {
    Ray oray = _worldToObj * ray;
    LocalGeo olocal;
    if (!_shape->intersect(oray, thit, &olocal)) return false;
    in->set_primitive(this);
    in->set_local(_objToWorld * olocal);
    return true;
}

bool GeometricPrimitive::intersectP(Ray& ray) {
    Ray oray = _worldToObj * ray;
    return _shape->intersectP(oray);
}

void GeometricPrimitive::getBRDF(LocalGeo& local, BRDF* brdf) {
    _mat->getBRDF(local, brdf);
}

BoundingBox& GeometricPrimitive::createBoundingBox() {
    return _shape->createBoundingBox();
}
/**END IMPLEMENTATION FOR GEOMETRIC PRIMITIVE. **/

/**BEGIN IMPLEMENTATION FOR AGGREGATE PRIMITIVE. **/
AggregatePrimitive::AggregatePrimitive() {
    //Default Constructor.
}

AggregatePrimitive::AggregatePrimitive(std::vector<Primitive*> list) {
    _list = list;
}

bool AggregatePrimitive::intersect(Ray& ray, float* thit, Intersection* in) {
    bool hit = false;
    for (int i = 0; i < _list.size(); i++) {
        if (_list.at(i)->intersect(ray, thit, in)) {
            ray.update_tmax(*thit);
            hit = true;
        }
    }
    return hit;
}

bool AggregatePrimitive::intersectP(Ray& ray) {
    for (int i = 0; i < _list.size(); i++) {
        if (_list.at(i)->intersectP(ray)) return true;
    }
    return false;
}

void AggregatePrimitive::getBRDF(LocalGeo& local, BRDF* brdf) {
    exit(1);
    //Should never enter here.
}

BoundingBox& AggregatePrimitive::createBoundingBox() {
    exit(1);
    //will implement later.
}

HBB::HBB() {
    //Default Constructor.
}

HBB::HBB(vector<Primitive*> list, int axis) {
    int len = list.size();
    if (len == 1) {
        _left = list.at(0);
        _right = NULL;
        _bbox = list.at(0)->createBoundingBox();
    } else if (len == 2) {
        _left = list.at(0);
        _right = list.at(1);
        _bbox = *new BoundingBox(list.at(0)->createBoundingBox(),
           list.at(1)->createBoundingBox());
    } else {
        BoundingBox listbox = makeListBox(list);
        float center = listbox.getCenter(axis);
        vector<Primitive*>* firstHalf = new vector<Primitive*>();
        vector<Primitive*>* secondHalf = new vector<Primitive*>();
        splitSpace(list, firstHalf, secondHalf, center, axis);
        _left = new HBB(*firstHalf, (axis + 1) % 3);
        _right = new HBB(*secondHalf, (axis + 1) % 3);
        _bbox = *new BoundingBox(_left->createBoundingBox(),
           _right->createBoundingBox());
    }
}

BoundingBox& HBB::makeListBox(vector<Primitive*> list) {
    BoundingBox* bb = new BoundingBox(list.at(0)->createBoundingBox(),
        list.at(0)->createBoundingBox());
    for (int i = 0; i < list.size(); i++) {
        bb = new BoundingBox(*bb, list.at(i)->createBoundingBox());
    }
    return *bb;
}

bool HBB::intersect(Ray& ray, float* thit, Intersection* in) {
    if (_bbox.intersectP(ray)) {
        Intersection* leftIn;
        Intersection* rightIn;
        float* leftT;
        float* rightT;
        bool leftHit = _left != NULL;
        if (leftHit) {
           leftHit = _left->intersect(ray, leftT, leftIn);
        }
        bool rightHit = _right != NULL;
        if (rightHit) {
            rightHit = _right->intersect(ray, rightT, rightIn);
        }
        if (leftHit && rightHit) {
            if (*leftT < *rightT) {
                thit = leftT;
                in = leftIn;
            } else {
                thit = rightT;
                in = rightIn;
            }
            return true;
        } else if (leftHit) {
            thit = leftT;
            in = leftIn;
            return true;
        } else if (rightHit) {
            thit = rightT;
            in = rightIn;
            return true;
        }
    }
    return false;
}

bool HBB::intersectP(Ray& ray) {
    if (_bbox.intersectP(ray)) {
        bool leftHit = _left != NULL;
        if (leftHit) {
           leftHit = _left->intersectP(ray);
        }
        bool rightHit = _right != NULL;
        if (rightHit) {
            rightHit = _right->intersectP(ray);
        }
        if (rightHit || leftHit) {
            return true;
        }
    }
    return false;
}

void HBB::getBRDF(LocalGeo& local, BRDF* brdf) {
    exit(1);
    //We should never enter here.
}

/** This method takes in LIST which is a list of primitives to be split.
    We split LIST on CENTER along AXIS. FIRSTHALF will contain all 
    objects less than CENTER. SECONDHALF will contain more centers
    more than CENTER. **/
void HBB::splitSpace(vector<Primitive*> list, vector<Primitive*>* firstHalf,
    vector<Primitive*>* secondHalf, float center, int axis) {
    for (int i = 0; i < list.size(); i++) {
        Primitive* primitive = list.at(i);
        float pCenter = primitive->createBoundingBox().getCenter(axis);
        if (pCenter >= center) {
            firstHalf->push_back(primitive);
        } else {
            secondHalf->push_back(primitive);
        }
    }
}

BoundingBox& HBB::createBoundingBox() {
    return _bbox;
}
/**END IMPLEMENTATION FOR AGGREGATE PRIMITIVE. **/

/**BEGIN IMPLEMENTATION FOR MATERIAL. **/
Material::Material(BRDF& brdf) {
    constantbrdf = brdf;
}

void Material::getBRDF(LocalGeo& local, BRDF* brdf) {
    brdf->setValue(constantbrdf);
}
/**END IMPLEMENTATION FOR MATERIAL. **/


/**BEGIN IMPLEMENTATION FOR SAMPLER. **/

Sampler::Sampler(int width, int height) {
    _width = width;
    _height = height;
    _currX = -1;
    _currY = 0;
}

Sampler::Sampler() {
    //empty constructor. 
}

bool Sampler::generateSample(Sample* sample) {
    _currX += 1;
    if (_currX == _width) {
        _currY += 1;
        _currX = 0;
    }
    if (_currY == _height) return false;
    sample->update_x(_currX);
    sample->update_y(_currY);
    return true;
}

/**END IMPLEMENTATION FOR SAMPLER. **/

/**BEGIN IMPLEMENTATION FOR CAMERA. **/
Camera::Camera() {
    //Empty constructor.
}

Camera::Camera(int width, int height, Point& pos, Point& looking,
    Vector& up, float fov) {
    _width = width;
    _height = height;
    _cameraPos = pos;
    Vector forward = looking - pos;
    _up = up;
    _up.normalize();
    _right = up.cross_product(forward);
    _right.normalize();
    _lookingAt = looking;
    float dif = abs(_cameraPos.get_z() - _lookingAt.get_z());
    t = dif * tan(fov / 2);
    b = -t;
    r = t * width / float(height);
    l = -r;
    LL = _lookingAt + _up * b + _right * l;
    LR = _lookingAt + _up * b + _right * r;
    UL = _lookingAt + _up * t + _right * l;
    UR = _lookingAt + _up * t + _right * r;
}

void Camera::generateRay(Sample& sample, Ray* ray) {
    float x = sample.get_x();
    float y = sample.get_y();
    float u = (x + .5) / _width;
    float v = (y + .5) / _height;
    Point p = (LL * v + UL * (1 - v)) * u  + (LR * v  + UR * (1 - v)) * (1 - u);
    Vector dir = p - _cameraPos;
    ray->update_origin(_cameraPos);
    ray->update_dir(dir);
    ray->update_tmin(0);
    ray->update_tmax(INFINITY);
}

Point& Camera::getCameraPos() {
    return _cameraPos;
}
/**END IMPLEMENTATION FOR CAMERA. **/

/**BEGIN IMPLEMENTATION FOR RAY TRACER. **/
RayTracer::RayTracer() {
    //default constructor
}

RayTracer::RayTracer(HBB& primitives, 
    std::vector<Light*> lights, Point& cameraPos) {
    _thit = INFINITY;
    _primitives = primitives;
    _lights = lights;
    _cameraPos = cameraPos;
}

void RayTracer::trace(Ray& ray, int depth, Color* color) {
    Color tempColor;
    BRDF brdf;
    Ray lray;
    Color lcolor;
    color->black();
    if (depth == 0) {
        color->black();
    } else {
        if (!_primitives.intersect(ray, &_thit, &_in)) {
            color->black();
        } else {
            _in.get_primitive()->getBRDF(_in.get_localGeo(), &brdf);
            for (int i = 0; i < _lights.size(); i++) {
                _lights.at(i)->generateLightRay(_in.get_localGeo(), &lray, &lcolor);
                if (!_primitives.intersectP(lray)) {
                    *color += shading(_in.get_localGeo(), brdf, lray, lcolor);
                }
            }

            if (brdf.isReflect()) {
                Ray reflectRay = createReflectRay(_in.get_localGeo(), ray);
                trace(reflectRay, depth-1, &tempColor);
                *color += brdf.get_kr() * tempColor;
            }
        }
    }
    _eye = _cameraPos;
}

Color& RayTracer::shading(LocalGeo& local, BRDF& brdf, 
    Ray& lray, Color& lcolor) {
    Color* shade = new Color();
    Normal norm = local.get_norm();
    Vector n;
    n.setValue(norm);
    Point p = local.get_pos();
    Vector l = lray.get_dir();
    Vector e = _eye - p;
    l.normalize();
    e.normalize();
    float lnCosine = l.dot_product(n);
    Vector k = (l * -1 + n * 2 * lnCosine);
    k.normalize();
    float keCosine = k.dot_product(e);
    *shade += brdf.get_kd() * max(lnCosine, 0);
    *shade += brdf.get_ka();
    *shade += brdf.get_ks() * pow(max(keCosine, 0), brdf.get_sp());
    *shade *= lcolor;
    return *shade;
}

Ray& RayTracer::createReflectRay(LocalGeo& local, Ray& ray) {
    float tmin = .01;
    float tmax = INFINITY;
    Point& point = local.get_pos();
    Normal norm = local.get_norm();
    Vector n;
    n.setValue(norm);
    Vector dir = ray.get_dir();
    dir = dir * -1;
    float cosine = dir.dot_product(n);
    Vector reflectDir = dir * -1 + n * 2 * cosine;
    Ray* reflectRay = new Ray(point, reflectDir, tmin, tmax);
    _eye = point;
    return *reflectRay;
}
/**END IMPLEMENTATION FOR RAY TRACER. **/

/**BEGIN IMPLEMENTATION FOR LIGHTS. **/
DirectionalLight::DirectionalLight() {
    //Default Constructor;
}

DirectionalLight::DirectionalLight(Vector& dir, Color& intensity) {
    _dir = dir;
    _intensity = intensity;
}

void DirectionalLight::generateLightRay(LocalGeo& local, Ray* ray, Color* lcolor) {
    ray->update_origin(local.get_pos());
    ray->update_dir(_dir * -1);
    ray->update_tmax(INFINITY);
    ray->update_tmin(.1);
    lcolor->setValue(_intensity);
}

PointLight::PointLight() {
    //Default Constructor;
}

PointLight::PointLight(Point& pos, Color& intensity) {
    _pos = pos;
    _intensity = intensity;
}

void PointLight::generateLightRay(LocalGeo& local, Ray* ray, Color* lcolor) {
    ray->update_origin(local.get_pos());
    ray->update_dir(_pos - local.get_pos());
    float maxt = _pos.distance(local.get_pos());
    ray->update_tmax(maxt);
    ray->update_tmin(0.1);
    lcolor->setValue(_intensity);
}
/**END IMPLEMENTATION FOR LIGHTS. **/

/**BEGIN IMPLEMENTATION FOR FILM. **/
Film::Film() {
    //Default constructor.
}

Film::Film(CImg<float>& img, bool display, const char* filename) {
    _img = img;
    _filename = filename;
    _display = display;
    if (display) {
        CImgDisplay main_display(img);
    }
}

void Film::commit(Sample& sample, Color& color) {
    int x = sample.get_x();
    int y = sample.get_y();
    const float pixel[] = {color.get_r(), color.get_g(), color.get_b()};
    _img.draw_point(x, y, pixel);
    if (_display) {
        CImgDisplay main_display(_img);
    }
}

void Film::writeImage() {
    LOG(INFO) << _filename;
    _img = scaleForJPG(_img);
    _img.save_jpeg(_filename);
}
/**END IMPLEMENTATION FOR FILM. **/

/**BEGIN IMPLEMENTATION FOR SCENE.  **/
Scene::Scene() {
    //Empty constructure
}

Scene::Scene(Sampler& sampler, Camera& camera,
    RayTracer& raytracer, Film& film, int depth) {
    _sampler = sampler;
    _camera = camera;
    _depth = depth;
    _raytracer = raytracer;
    _film = film;
}

void Scene::render() {
    while(_sampler.generateSample(&_sample)) {
        _camera.generateRay(_sample, &_ray);
        _raytracer.trace(_ray, _depth, &_color);
        _film.commit(_sample, _color);
    }
    _film.writeImage();
}


/**END IMPLEMENTATION FOR SCENE. **/

