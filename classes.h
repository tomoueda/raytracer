#include <vector>
#include <stdlib.h>

class Vector;
class Normal;
class LocalGeo;
class Shape;
class Material;
class Intersection;

class Vector {
    friend class Point;
    friend class Color;
    float _x, _y, _z;
public:
    Vector();
    Vector(float x, float y, float z);
    Vector & operator+(Vector &v);
    void operator+=(Vector &v);
    Vector & operator-(Vector &v);
    void operator-=(Vector &v);
    Vector & operator*(float scalar);
    void operator*=(float scalar);
    Vector & operator/(float scalar);
    void operator/=(float scalar);
    bool operator==(const Vector &v) const;
    bool operator!=(const Vector &v) const;
    Normal & normalize();
    void setValue(float x, float y, float z);
    void debug();
};

class Normal {
    float _x, _y, _z;
public:
    Normal();
    Normal(float x, float y, float z);
    void normalize();
    Normal & operator+(Normal &n);
    Normal & operator-(Normal &n);
    void operator+=(Normal &n);
    void operator-=(Normal &n);
    bool operator==(const Normal& n) const;
    bool operator!=(const Normal& n) const;
    void debug();
};

class Point {
    float _x, _y, _z;
public:
    Point();
    Point(float x, float y, float z);
    Point & operator+(Vector &v);
    Point & operator-(Vector &v);
    void operator+=(Vector &v);
    void operator-=(Vector &v);
    Vector & operator - (Point &p);
    bool operator==(const Point& p) const;
    bool operator!=(const Point& p) const;
    void debug();
};

class Ray {
    Point _pos;
    Vector _dir;
    float _t_min, _t_max;
public:
    Ray(Point &pos, Vector &dir, float t_min, float t_max);
    float get_t_min();
    float get_t_max();
};

class Matrix {
    float mat[4][4];
};

class Rotation: public Matrix {
public:
    Rotation(Point &about, float degree);
};

class Translation: public Matrix {
public:
    Translation(float x, float y, float z);
};

class Scaling: public Matrix {
public:
    Scaling(float x, float y, float z);
};

class Transformation {
    Matrix *m;
public:
    Point & operator*(Point &p);
    Vector & operator*(Vector &v);
    Normal & operator*(Normal &n);
    Ray & operator*(Ray &r);
    LocalGeo & operator*(LocalGeo *l);
};

class Color {
    float _r, _g, _b;
public:
    Color();
    Color(float r, float g, float b);
    Color &operator+(Color &c);
    Color &operator-(Color &c);
    Color &operator*(Color &c);
    Color &operator/(Color &c);
    void operator+=(Color &c);
    void operator-=(Color &c);
    void operator*=(Color &c);
    void operator/=(Color &c);
    Color &operator*(float s);
    Color &operator/(float s);
    void operator*=(float s);
    void operator/=(float s);
    void setValue(float r, float g, float b);
    bool operator==(const Color& c) const;
    bool operator!=(const Color& c) const;
    void debug();
    void convert(Vector &v);
};

class BRDF {
    Color _kd, _ks, _ka, _kr;
public:
    BRDF(Color &kd, Color &ks, Color &ka, Color &kr);
    Color& get_kd();
    Color& get_ks();
    Color& get_ka();
    Color& get_kr();
};

class Sample {
    float _x, _y;
public:
    Sample(float x, float y);
    float get_x();
    float get_y();
};

class LocalGeo {
    Point _pos;
    Normal _normal;
public:
    LocalGeo(Point &pos, Point &normal);
    Point& get_pos();
    Normal& get_norm();
};

/***** Other Classes *****/

class Primitive {
public:
    virtual bool intersect(Ray &ray, float *thit, Intersection *in);
    virtual bool intersectP(Ray &ray);
    virtual void getBRDF(LocalGeo &local, BRDF *brdf);
};

class Intersection {
    LocalGeo _localGeo;
    Primitive* _primitive;
public:
    Intersection(LocalGeo &localGeo, Primitive *primitive);
    LocalGeo& 
};

class GeometricPrimitive: public Primitive {
    Transformation objToWorld, worldToObj;
    Shape *shape;
    Material *mat;
public:
    bool intersect(Ray &ray, float *thit, Intersection *in);
    bool intersectP(Ray &ray);
    void getBRDF(LocalGeo &local, BRDF *brdf);
};

class AggregatePrimitive: public Primitive {
public:
    AggregatePrimitive(std::vector<Primitive*> list);
    bool intersect(Ray &ray, float *thit, Intersection *in);
    bool intersectP(Ray &ray);
    void getBRDF(LocalGeo &local, BRDF *brdf);
};

class Material {
    BRDF const brdf;
public:
    void getBRDF(LocalGeo &local, BRDF*brdf);
};

class Sampler {
public:
    bool getSample(Sample* sample);
};

class Camera {
public:
    void generateRay(Sample &sample, Ray *ray);
};

class RayTracer {
public:
    void trace(Ray &ray, int depth, Color *color);
};

class Light {
public:
    virtual void generateLightRay(LocalGeo &local, Ray *ray, Color *lcolor);
};

class DirectionalLight: public Light {
public:
    void generateLightRay(LocalGeo &local, Ray *ray, Color *lcolor);
};

class PointLIght: public Light {
public:
    void generateLightRay(LocalGeo &local, Ray *ray, Color *lcolor);
};

class Film {
public:
    void commit(Sample &sample, Color &color);
    void writeImage();
};

class Scene {
public:
    void render();
};
