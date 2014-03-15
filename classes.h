#include <vector>
#include <stdlib.h>
#include <CImg.h>
#undef Success
#include <Eigen/Dense>

using namespace cimg_library;
using namespace Eigen;
using namespace std;

class Vector;
class Normal;
class LocalGeo;
class Shape;
class Material;
class Intersection;
class Light;
class BoundingBox;
class Primitive;



class Vector {
    friend class Point;
    friend class Color;
    friend class Normal;
    float _x, _y, _z;
public:
    Vector();
    Vector(float x, float y, float z);
    Vector & operator+(Vector &v);
    Vector & operator+(Normal &v);
    void operator+=(Vector &v);
    Vector & operator-(Vector &v);
    void operator-=(Vector &v);
    Vector & operator*(float scalar);
    void operator*=(float scalar);
    Vector & operator/(float scalar);
    void operator/=(float scalar);
    bool operator==(const Vector &v) const;
    bool operator!=(const Vector &v) const;
    Vector& cross_product(Vector& v); 
    Normal& get_normal();
    void normalize();
    void setValue(float x, float y, float z);
    void setValue(Normal& n);
    float dot_product(Vector& v);
    float dot_product(Normal& v);
    float get_x();
    float get_y();
    float get_z();
    void debug();
};

class Normal {
    friend class Point;
    friend class Vector;
    float _x, _y, _z;
public:
    Normal();
    Normal(float x, float y, float z);
    void normalize();
    Normal & operator+(Normal &n);
    Normal & operator-(Normal &n);
    void operator+=(Normal &n);
    void operator-=(Normal &n);
    Normal & operator*(float s);
    bool operator==(const Normal& n) const;
    bool operator!=(const Normal& n) const;
    Vector& convertToVector();
    float get_x();
    float get_y();
    float get_z();
    float dot_product(Normal& n);
    float dot_product(Vector& v);
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
    Point& operator*(float f);
    Point& operator+(Point& p);
    bool operator==(const Point& p) const;
    bool operator!=(const Point& p) const;
    float dot_product(Normal& n);
    float dot_product(Point& p);
    float get_x();
    float get_y();
    float get_z();
    void debug();
    float distance(Point& p);
};

class Ray {
    Point _pos;
    Vector _dir;
    float _t_min, _t_max;
public:
    Ray();
    Ray(Point& pos, Vector& dir, float t_min, float t_max);
    float get_t_min();
    float get_t_max();
    Point& get_pos();
    Vector& get_dir();
    Point& get_pos_with_t(float t);
    void update_tmax(float t);
    void update_origin(Point& o);
    void update_dir(Vector& dir);
    void update_tmin(float t);
};

class TransMatrix {
    Matrix4f _mat;
public:
    TransMatrix();
    TransMatrix(Matrix4f mat);
    void add_rotation(Point& about, float radian);
    void add_translation(float x, float y, float z);
    void add_scaling(float x, float y, float z);
    TransMatrix& invt();
    TransMatrix& inv();
    Vector4f operator*(Vector4f v);
};

class Transformation {
    TransMatrix _m;
    TransMatrix _minvt;
public:
    Transformation();
    Transformation(TransMatrix& m);
    Point& operator*(Point &p);
    Vector& operator*(Vector &v);
    Normal& operator*(Normal &n);
    Ray& operator*(Ray &r);
    LocalGeo& operator*(LocalGeo& l);
    Transformation& inv();
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
    void setValue(Color& c);
    void black();
    void white();
    bool operator==(const Color& c) const;
    bool operator!=(const Color& c) const;
    void debug();
    void convert(Vector &v);
    float get_r();
    float get_g();
    float get_b();
};

class BRDF {
    Color _kd, _ks, _ka, _kr;
    float _sp;
public:
    BRDF();
    BRDF(Color &kd, Color &ks, Color &ka, Color &kr, float sp);
    Color& get_kd();
    Color& get_ks();
    Color& get_ka();
    Color& get_kr();
    float get_sp();
    void setValue(BRDF& brdf);
    bool isReflect();
};

class Sample {
    float _x, _y;
public:
    Sample();
    Sample(float x, float y);
    float get_x();
    float get_y();
    void update_x(float x);
    void update_y(float y);
};

class LocalGeo {
    Point _pos;
    Normal _normal;
public:
    LocalGeo();
    LocalGeo(Point& pos, Normal& normal);
    Point& get_pos();
    Normal& get_norm();
    void set_pos(Point& p);
    void set_norm(Normal& n);
};

/***** Other Classes *****/
class Vertex {
    Normal _norm;
    Point _pos;
public:
    Vertex();
    Vertex(Point& pos, Normal& norm);
    Normal& get_norm();
    Point& get_point();
    void set_norm(Normal& norm);
    void set_pos(Point& pos);
};

class Shape {
public:
    virtual bool intersect(Ray& ray, float* thit, LocalGeo* local) = 0;
    virtual bool intersectP(Ray& ray) = 0;
    virtual BoundingBox& createBoundingBox() = 0;
};

class Sphere: public Shape {
    Point _center;
    float _radius;
public:
    Sphere(Point& center, float radius);
    bool intersect(Ray& ray, float* thit, LocalGeo* local);
    bool intersectP(Ray& ray);
    Point& getCenter();
    float getRadius();
    BoundingBox& createBoundingBox();
};

class Triangle: public Shape {
    Point _v1, _v2, _v3;
public:
    Triangle(Point& v1, Point& v2, Point& v3);
    bool intersect(Ray& ray, float* thit, LocalGeo* local);
    bool intersectP(Ray& ray);
    Point& getv1();
    Point& getv2();
    Point& getv3();
    BoundingBox& createBoundingBox();
};

class VertexTriangle: public Shape {
    Vertex _v1, _v2, _v3;
public:
    VertexTriangle(Vertex& v1, Vertex& v2, Vertex& v3);
    bool intersect(Ray& ray, float* thit, LocalGeo* local);
    bool intersectP(Ray& ray);
    Vertex& getv1();
    Vertex& getv2();
    Vertex& getv3();
    BoundingBox& createBoundingBox();
};

class BoundingBox: public Shape {
    float _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;
    Point _max, _min;
public:
    BoundingBox();
    BoundingBox(BoundingBox& bb1, BoundingBox& bb2);
    BoundingBox(Point& max, Point& min);
    BoundingBox(float xmax, float ymax, float zmax,
        float xmin, float ymin, float zmin);
    bool intersect(Ray& ray, float* thit, LocalGeo* local);
    bool intersectP(Ray& ray);
    BoundingBox& createBoundingBox();
    float getCenter(int axis);
    void extend(BoundingBox& bb);
    void debug();
    Point& getMax();
    Point& getMin();
};

class Primitive {
public:
    virtual bool intersect(Ray &ray, float *thit, Intersection *in) = 0;
    virtual bool intersectP(Ray &ray) = 0;
    virtual void getBRDF(LocalGeo &local, BRDF *brdf) = 0;
    virtual BoundingBox& createBoundingBox() = 0;
    virtual void debug() = 0;
};

class Intersection {
    LocalGeo _localGeo;
    Primitive* _primitive;
public:
    Intersection();
    Intersection(LocalGeo &localGeo, Primitive *primitive);
    LocalGeo& get_localGeo();
    Primitive* get_primitive();
    void set_primitive(Primitive* primitive);
    void set_local(LocalGeo& localGeo);
};

class GeometricPrimitive: public Primitive {
    Transformation _objToWorld, _worldToObj;
    Shape* _shape;
    Material* _mat;
public:
    GeometricPrimitive();
    GeometricPrimitive(Transformation& trans, Shape* shape, Material* mat);
    bool intersect(Ray &ray, float *thit, Intersection *in);
    bool intersectP(Ray &ray);
    void getBRDF(LocalGeo &local, BRDF *brdf);
    BoundingBox& createBoundingBox();
    void debug();
};

class AggregatePrimitive: public Primitive {
    std::vector<Primitive*> _list;
public:
    AggregatePrimitive();
    AggregatePrimitive(std::vector<Primitive*> list);
    bool intersect(Ray &ray, float *thit, Intersection *in);
    bool intersectP(Ray &ray);
    void getBRDF(LocalGeo &local, BRDF *brdf);
    BoundingBox& createBoundingBox();
    void debug();
};

class HBB : public Primitive {
    Primitive* _left;
    Primitive* _right;
    BoundingBox _bbox;
public:
    HBB();
    HBB(vector<Primitive*> list, int axis);
    bool intersect(Ray& ray, float* thit, Intersection* in);
    bool intersectP(Ray& ray);
    void getBRDF(LocalGeo& local, BRDF* brdf);
    void splitSpace(vector<Primitive*> list,
        vector<Primitive*>* firstHalf, vector<Primitive*>* secondHalf,
        float center, int axis);
    BoundingBox& createBoundingBox();
    BoundingBox& makeListBox(vector<Primitive*> list);
    void debug();
};

class Material {
    BRDF constantbrdf;
public:
    Material(BRDF& brdf);
    void getBRDF(LocalGeo &local, BRDF*brdf);
};

class Sampler {
    int _width;
    int _height;
    float _currX;
    float _currY;
public:
    Sampler();
    Sampler(int width, int height);
    bool generateSample(Sample* sample);
};

class Camera {
    int _width;
    int _height;
    Point _cameraPos, _lookingAt, LL, UL, UR, LR;
    Vector _up;
    Vector _right;
    float t, b, r, l;
public:
    Camera();
    Camera(int width, int height, Point& pos, Point& looking,
        Vector& up, float fov);
    void generateRay(Sample& sample, Ray *ray);
    Point& getCameraPos();
};

class RayTracer {
    float _thit;
    Intersection _in;
    HBB _primitives;
    AggregatePrimitive _list;
    Ray _lray;
    std::vector<Light*> _lights;
    Color _lcolor;
    Point _cameraPos;
    Point _eye;
public:
    RayTracer();
    RayTracer(HBB& primitives,
        std::vector<Light*> lights, Point& cameraPos);
    RayTracer(AggregatePrimitive& primitives,
        std::vector<Light*> lights, Point& cameraPos);
    void trace(Ray &ray, int depth, Color *color);
    void traceUsingList(Ray &ray, int depth, Color *color);
    Color& shading(LocalGeo& local, BRDF& brdf, Ray& lray, Color& lcolor);
    Ray& createReflectRay(LocalGeo& geo, Ray& ray);
};

class Light {
protected:
    Color _intensity;
public:
    virtual void generateLightRay(LocalGeo &local, Ray *ray, Color *lcolor)=0;
};

class DirectionalLight: public Light {
    Vector _dir;
public:
    DirectionalLight();
    DirectionalLight(Vector& dir, Color& intensity);
    void generateLightRay(LocalGeo &local, Ray *ray, Color *lcolor);
};

class PointLight: public Light{
    Point _pos;
public:
    PointLight();
    PointLight(Point& pos, Color& intensity);
    void generateLightRay(LocalGeo &local, Ray *ray, Color *lcolor);
};

class Film {
    CImg<float> _img;
    const char* _filename;
    bool _display;
public:
    Film();
    Film(CImg<float>& img, bool display, const char* filename);
    void commit(Sample &sample, Color &color);
    void writeImage();
};

class Scene {
    Sampler _sampler;
    Camera _camera;
    RayTracer _raytracer;
    Film _film;
    Sample _sample;
    Ray _ray;
    Color _color;
    int _depth;
public:
    Scene();
    Scene(Sampler& sampler, Camera& camera,
       RayTracer& raytracer, Film& film, int depth);
    void render();
    void renderUsingList();
};
