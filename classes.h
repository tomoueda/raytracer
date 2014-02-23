#include <vector>

class Vector;
class Normal;
class LocalGeo;
class Shape;
class Material;
class Intersection;



class Vector {
    float _x, _y, _z;
public:
    Vector(float x, float y, float z);
    Vector & operator+(Vector &v);
    Vector & operator-(Vector &v);
    Vector & operator*(float scalar);
    Vector & operator/(float scalar);
    Normal & normalize();
};

class Normal {
    float x, y, z;
public:
    Normal(float x, float y, float z);
    Normal & operator+(Normal &n);
    Normal & operator-(Normal &n);
};

class Point {
    float x, y, z;
public:
    Point(float x, float y, float z);
    Point & operator+(Vector &v);
    Point & operator-(Vector &v);
    Vector & operator - (Point &p);
};

class Ray {
    Point pos;
    Vector dir;
    float t_min, t_max;
public:
    Ray(Point &pos, Vector &dir, float t_min, float t_max);
};

class Matrix {
    float mat[4][4];
    void invert();
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
    float r, g, b;
public:
    Color &operator+(Color &c);
    Color &operator-(Color &c);
    Color &operator*(float s);
    Color &operator/(float s);
};

class BRDF {
    Color *kd;
    Color *ks;
    Color *ka;
    Color *kr;
public:
    BRDF(Color &kd, Color &ks, Color &ka, Color &kr);
};

class Sample {
    float x, y;
public:
    Sample(float x, float y);
};

class LocalGeo {
    Point *pos;
    Normal *normal;
public:
    LocalGeo(Point &pos, Point &normal);
};

/***** Other Classes *****/

class Primitive {
public:
    virtual bool intersect(Ray &ray, float *thit, Intersection *in);
    virtual bool intersectP(Ray &ray);
    virtual void getBRDF(LocalGeo &local, BRDF *brdf);
};

class Intersection {
    LocalGeo localGeo;
    Primitive *primitive;
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
