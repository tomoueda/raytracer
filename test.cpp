#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include "CImg.h"
#include "classes.h"
#include "as2.h"

#define pi atan(1)*4

using namespace std;

//Tests for the Vector Class
TEST(classes, vector) {
	Vector a(1.1, 2.2, 3.3);
	Vector b(4, 5, -9);
	Vector d(1, 1, 1);
	d = a + b;
	Vector c(5.1, 7.2, -5.7);
	Vector e(1.1, 2.2, 3.3);
	EXPECT_EQ(d, c) << "testing '+' operator";
	EXPECT_EQ(a, e) << "Making sure that a is untouched. ";
	a += b;
	EXPECT_EQ(d, a) << "Testing out '+=' operator";
	a.setValue(1, 2, 3);
	a = a - b;
	c.setValue(-3, -3, 12);
	EXPECT_EQ(a, c) << "testing '-' operator ";
	c -= b;
	a.setValue(-7, -8, 21);
	EXPECT_EQ(a, c) << "testing '-=' operator";
	a.setValue(1, 2, 3);
	float scale = 5;
	a = a * scale;
	c.setValue(5, 10, 15);
	EXPECT_EQ(a, c) << "testing '*' operator";
	a.setValue(1, 2, 3);
	a *= 5;
	EXPECT_EQ(a, c) << "testing '*=' operator ";
	a.setValue(1, 2, 3);
	a = a / scale;
	c.setValue(1.0/5, 2.0/5, 3.0/5);
	EXPECT_EQ(a, c) << "testing '/' operator";
	a.setValue(1, 2, 3);
	a /= scale;
	EXPECT_EQ(a, c) << "testing '/=' operator";
	c.setValue(10, 10, 10);
	EXPECT_NE(a, c) << "testing the '!=' operator";
	a.setValue(10, 10, 10);
	c.setValue(1, 2, 3);
	float test = a.dot_product(c);
	EXPECT_EQ(test, 60) << "Testing dot product.";
}

TEST(classes, normal) {
	Vector a(.4, 0, 0);
	Normal q(.1230921, 0, 0);
	Normal b(1, 0, 0);
	Normal c = a.get_normal();
	EXPECT_EQ(c, b) << "testing normalize";
	Normal d(.000001, 0, 0);
	Normal e(0, 0, 0);
	Normal f(5, 5, 5);
	Normal g(.57735, .57735, .57735);
	EXPECT_EQ(d, e) << "testing constructor normalize roundoff";
	EXPECT_EQ(f, g) << "testing constructor normalize";
	d = d+q;
	EXPECT_EQ(b, d) << "testing '+' operator";
	Normal j(.0000001, 0, 0);
	j += q;
	EXPECT_EQ(j, q) << "testing '+=' operator";
	b = b+ f;
	Normal k(.888074, .325058, .325058);
	EXPECT_EQ(k, b) << "case 5";
}

TEST(classes, point) {
	Point a(1, 2, 3);
	Point b(10, 10, 10);
	Vector z(-9, -8, -7);
	Vector y = a - b;
	EXPECT_EQ(y, z) << "Point - Point = Vector";
	a = a + z;
	Point c(-8, -6, -4);
	EXPECT_EQ(a, c) << "testing '+' operator";
	c = c - z;
	Point q(1, 2, 3);
	EXPECT_EQ(c, q) << "testing '-' operator";

}

//Tests for the Vector Class
TEST(classes, color) {
	Color a(1.1, 2.2, 3.3);
	Color b(4, 5, -9);
	Color d(1, 1, 1);
	d = a + b;
	Color c(5.1, 7.2, -5.7);
	Color e(1.1, 2.2, 3.3);
	EXPECT_EQ(d, c) << "testing '+' operator";
	EXPECT_EQ(a, e) << "Making sure that a is untouched. ";
	a += b;
	EXPECT_EQ(d, a) << "Testing out '+=' operator";
	a.setValue(1, 2, 3);
	a = a - b;
	c.setValue(-3, -3, 12);
	EXPECT_EQ(a, c) << "testing '-' operator ";
	c -= b;
	a.setValue(-7, -8, 21);
	EXPECT_EQ(a, c) << "testing '-=' operator";
	a.setValue(1, 2, 3);
	float scale = 5;
	a = a * scale;
	c.setValue(5, 10, 15);
	EXPECT_EQ(a, c) << "testing '*' operator";
	a.setValue(1, 2, 3);
	a *= 5;
	EXPECT_EQ(a, c) << "testing '*=' operator ";
	a.setValue(1, 2, 3);
	a = a / scale;
	c.setValue(1.0/5, 2.0/5, 3.0/5);
	EXPECT_EQ(a, c) << "testing '/' operator";
	a.setValue(1, 2, 3);
	a /= scale;
	EXPECT_EQ(a, c) << "testing '/=' operator";
	c.setValue(10, 10, 10);
	EXPECT_NE(a, c) << "testing the '!=' operator";
	a.setValue(1, 2, 3);
	a *= a;
	Color tester(1, 4, 9);
	EXPECT_EQ(a, tester) << "testing the 'vector mult'";
	a /= a;
	tester.setValue(1, 1, 1);
	EXPECT_EQ(a, tester) << "testing the 'vector div'";
	Vector vec(12.34, 12.2, 120.9);
	a.convert(vec);
	tester.setValue(12.34, 12.2, 120.9);
	EXPECT_EQ(a, tester) << "testing the vector conversion";
}

TEST(classes, brdf) {
	Color kd(.8, .5, .3);
	Color ks(.9, .9, .8);
	Color ka(.1, .2, .4);
	Color kr(.8, .5, .9);
	BRDF brdf(kd, ks, ka, kr, 19);
	EXPECT_EQ(kd, brdf.get_kd()) << "test kd equality";
	EXPECT_EQ(ks, brdf.get_ks()) << "test ks equality";
	EXPECT_EQ(ka, brdf.get_ka()) << "test ka equality";
	EXPECT_EQ(kr, brdf.get_kr()) << "test kr equality";
}

TEST(classes, sample) {
	Sample sample(34.566, 9844);
	float x = 34.566;
	float y = 9844;
	EXPECT_EQ(sample.get_x(), x) << "testing get x";
	EXPECT_EQ(sample.get_y(), y) << "testing get y";
}

TEST(classes, sphere_trivial) {
	Point center(0, 0, 0);
	Point ray_origin(0, 0, 5);
	Vector ray_dir(0, 0, -1);
	float tmin = 0;
	float tmax = INFINITY;
	float radius = 1;
	float* t;
	LocalGeo* localGeo = new LocalGeo();
	Sphere sphere(center, radius);
	Ray ray(ray_origin, ray_dir, tmin, tmax);
	bool intersected = sphere.intersect(ray, t, localGeo);
	Point localGeo_p = localGeo->get_pos();
	Normal localGeo_n = localGeo->get_norm();
	Point testp(0, 0, 1);
	Normal testn(0, 0, 1);
	float test_t = 4;
	/** This is a trivial test, for sanity check. **/
	EXPECT_EQ(intersected, true) << "Testing if it intersected.";
	EXPECT_EQ(localGeo_p, testp) << "Testing Local Geo position update.";
	EXPECT_EQ(localGeo_n, testn) << "Testing Local Geo norm update.";
	EXPECT_EQ(test_t, *t) << "Testing the t update.";
}

TEST(classes, sphere_nontrivial) {
	Point center(-12, 11, 18);
	Point ray_origin(12, -11, -18);
	Vector ray_dir(1, 1, 1);
	Vector ray_dir2(-24, 22, 36);
	Vector ray_dir3(0, 0, 0);
	float tmin = 0;
	float inftmin = INFINITY;
	float tmax = INFINITY;
	float smalltmax = 0;
	float radius = 5;
	float* t;
	LocalGeo* localGeo = new LocalGeo();
	Sphere sphere(center, radius);
	Ray ray(ray_origin, ray_dir, tmin, tmax);
	Ray rayhit(ray_origin, ray_dir2, tmin, tmax);
	Ray bigtmin(ray_origin, ray_dir2, inftmin, tmax);
	Ray smallesttmax(ray_origin, ray_dir2, tmin, smalltmax);
	Ray nonmoving(ray_origin, ray_dir3, tmin, tmax);
	bool intersected = sphere.intersect(ray, t, localGeo);
	bool intersected2 = sphere.intersectP(rayhit);
	bool intersected3 = sphere.intersectP(bigtmin);
	bool intersected4 = sphere.intersectP(smallesttmax);
	bool intersected5 = sphere.intersectP(nonmoving);
	EXPECT_NE(intersected, true) << "Testing if a ray misses";
	EXPECT_EQ(true, intersected2) << "Testing intersectP";
	EXPECT_EQ(false, intersected3) << "Testing tmin check.";
	EXPECT_EQ(false, intersected4) << "Testing tmax check.";
	EXPECT_EQ(false, intersected5) << "Testing nonmoving ray";
}

TEST(classes, triangle) {
	Point v1(1, 2, 0);
	Point v2(-1, 0, 0);
	Point v3(1, -2, 0);
	Triangle tri(v1, v2, v3);
	Point ray_origin(0, 0, 5);
	Vector ray_dir(0, 0, -1);
	Vector ray_dir2(0, 0, 0);
	float tmin = 0;
	float tmax = INFINITY;
	float bigmin = INFINITY;
	float smallmax = 3;
	float* t;
	LocalGeo* localGeo = new LocalGeo();
	Ray ray(ray_origin, ray_dir, tmin, tmax);
	Ray largemin(ray_origin, ray_dir, bigmin, tmax);
	Ray smallmaxi(ray_origin, ray_dir, tmin, smallmax);
	Ray nonmoving(ray_origin, ray_dir2, tmin, tmax);
	bool intersected = tri.intersect(ray, t, localGeo);
	bool intersected2 = tri.intersectP(largemin);
	bool intersected3 = tri.intersectP(smallmaxi);
	bool intersected4 = tri.intersectP(nonmoving);
	Point inter(0, 0, 0);
    Normal norm	= ((v2 - v1).cross_product(v3 - v1)).get_normal();
    float test = 5;
    EXPECT_EQ(true, intersected) << "making sure it intersects";
    EXPECT_EQ(test, *t) << "Testing the t update";
    EXPECT_EQ(inter, localGeo->get_pos()) << "Testing the pos update.";
    EXPECT_EQ(norm, localGeo->get_norm()) << "Testing the norm update.";
    EXPECT_EQ(false, intersected2) << "Testing intersectP and large min";
    EXPECT_EQ(false, intersected3) << "Testing intersectP and small max. ";
    EXPECT_EQ(false, intersected4) << "Testing nonmoving ray";
}

TEST(classes, transformation) {
	TransMatrix test;
	Point a(1, 2, 3);
	Vector b(2, 3, 4);
	Normal c(1, 0, 0);
	Ray d(a, b, 0, 10);
	LocalGeo e(a, c);
	test.add_translation(4, 5, 9);
	Transformation t(test);
	a = t * a;
	b = t * b;
	c = t * c;
	d = t * d;
	e = t * e;
	Point f(5, 7, 12);
	Vector g(2, 3, 4);
	Ray i(f, g, 0, 10);
	LocalGeo j(f, c);
	EXPECT_EQ(f, a) << "Testing point trans";
	EXPECT_EQ(g, b) << "Testing Vector trans";
	EXPECT_EQ(i.get_pos(), d.get_pos()) << "Testing the ray trans: pos";
	EXPECT_EQ(j.get_pos(), e.get_pos()) << "Testing the geolocal trans";
	Point z(1, 0, 0);
	Point cen(0, -1, 0);
	TransMatrix test2;
	test2.add_rotation(cen, pi/2);
	Transformation t2(test2);
	z = t2 * z;
	Point expect(0, 0, 1);
	EXPECT_EQ(expect, z) << "Testing rotation";
	TransMatrix test3;
	test3.add_rotation(cen, pi/2);
	test3.add_scaling(2, 1, 1);
	Transformation t3(test3);
	Point z2(1, 0, 0);
	z2 = t3 * z2;
	Point expect2(0, 0, 2);
	z2.debug();
	EXPECT_EQ(expect2, z2) << "Testing multiple transformation.";
	Transformation test3inv = t3.inv();
	z2 = test3inv * z2;
	Point expect3(1, 0, 0);
	EXPECT_EQ(expect3, z2) << "Testing inverse transformation.";
}

TEST(class, GeometricPrimitive) {
	TransMatrix transM;
	transM.add_translation(1, 2, 3);
	Transformation trans(transM);
	Point camera(1, 2, 10);
	Vector dir(0, 0, -50);
	Ray ray(camera, dir, 0, INFINITY);
	Point center(0, 0, 0);
	float radius = 1;
	Sphere* sphere = new Sphere(center, radius);
	Color kd(.9, .9, .4);
	Color ks(.3, .4, .5);
	Color ka(.1, .2, .4);
	Color kr(.4, .4, .5);
	float sp = 16;
	BRDF brdf(kd, ks, ka, kr, sp);
	Material* material = new Material(brdf);
	GeometricPrimitive prim(trans, sphere, material);
	Intersection* intersect = new Intersection();
	float* thit;
	bool intersected =  prim.intersect(ray, thit, intersect);
	BRDF* brdftester = new BRDF();
	prim.getBRDF(intersect->get_localGeo(), brdftester);
	float test = 6;
	LocalGeo tester = intersect->get_localGeo();
	Point inter(1, 2, 4);
	Normal norm(0, 0, 1);
	EXPECT_EQ(inter, tester.get_pos()) << "Testing the right intersection position.";
	EXPECT_EQ(norm, tester.get_norm()) << "Testing the right intersection normal.";
	EXPECT_EQ(brdftester->get_kd(), brdf.get_kd()) << "testing brdf copying";
	EXPECT_EQ(brdftester->get_ks(), brdf.get_ks()) << "testing brdf copying";
	EXPECT_EQ(brdftester->get_ka(), brdf.get_ka()) << "testing brdf copying";
	EXPECT_EQ(brdftester->get_kr(), brdf.get_kr()) << "testing brdf copying";
	EXPECT_EQ(brdftester->get_sp(), brdf.get_sp()) << "testing brdf copying";
	EXPECT_EQ(test, *thit) << "testing right value.";
	EXPECT_EQ(true, intersected) << "Testing if it even hit.";
}

TEST(class, AggregatePrimitive) {
	TransMatrix transM;
	transM.add_translation(1, 2, 3);
	Transformation trans(transM);
	Point camera(1, 2, 10);
	Vector dir(0, 0, -50);
	Ray ray(camera, dir, 0, INFINITY);
	Point center(0, 0, 0);
	float radius = 1;
	Point center2(0, 0, -1000);
	Sphere* sphere = new Sphere(center, radius);
	Sphere* sphere2 = new Sphere(center2, radius);
	Color kd(.9, .9, .4);
	Color ks(.3, .4, .5);
	Color ka(.1, .2, .4);
	Color kr(.4, .4, .5);
	float sp = 16;
	BRDF brdf(kd, ks, ka, kr, sp);

	Color kd2(100000, .9, .4);
	Color ks2(.3, .4, .5);
	Color ka2(.1, .2, 400000);
	Color kr2(.4, 0, .5);
	float sp2 = 16000;
	BRDF brdf2(kd2, ks2, ka2, kr2, sp2);

	Material* material = new Material(brdf);
	Material* material2 = new Material(brdf2);
	GeometricPrimitive* prim = new GeometricPrimitive(trans, sphere, material);
	GeometricPrimitive* prim2 = new GeometricPrimitive(trans, sphere2, material2);
	std::vector<Primitive*> list;
	list.push_back(prim);
	list.push_back(prim2);
	AggregatePrimitive aggs(list);
	LOG(INFO) << list.size();
	Intersection* intersect = new Intersection();
	float* thit;
	bool intersected =  aggs.intersect(ray, thit, intersect);
	bool intersectedP = aggs.intersectP(ray);
	BRDF* brdftester = new BRDF();
	intersect->get_primitive()->getBRDF(intersect->get_localGeo(), brdftester);
	float test = 6;
	LocalGeo tester = intersect->get_localGeo();
	Point inter(1, 2, 4);
	Normal norm(0, 0, 1);
	EXPECT_EQ(inter, tester.get_pos()) << "Testing the right intersection position.";
	EXPECT_EQ(norm, tester.get_norm()) << "Testing the right intersection normal.";
	EXPECT_EQ(brdftester->get_kd(), brdf.get_kd()) << "testing brdf copying";
	EXPECT_EQ(brdftester->get_ks(), brdf.get_ks()) << "testing brdf copying";
	EXPECT_EQ(brdftester->get_ka(), brdf.get_ka()) << "testing brdf copying";
	EXPECT_EQ(brdftester->get_kr(), brdf.get_kr()) << "testing brdf copying";
	EXPECT_EQ(brdftester->get_sp(), brdf.get_sp()) << "testing brdf copying";
	EXPECT_EQ(test, *thit) << "testing right value.";
	EXPECT_EQ(true, intersected) << "Testing if it even hit.";
	EXPECT_EQ(true, intersectedP) << "Testing if the p version works";
}

TEST(class, lights) {
	Point pos(0.0, 0.0, 0.0);
	Normal norm(0, 1, 0);
	LocalGeo geo(pos, norm);
	Ray* lray = new Ray();
	Color* lcolor = new Color();
	Color color(3, .4, .5);
	Point p(2.0, 2.0, 0.0);
	PointLight point(p, color);
	point.generateLightRay(geo, lray, lcolor);
	Vector dir(2, 2, 0);
	dir.normalize();
	EXPECT_EQ(pos, lray->get_pos()) << "Testing lray position";
	EXPECT_EQ(dir, lray->get_dir()) << "Getting lray direction";
	EXPECT_EQ(color, *lcolor) << "Make sure the color translated";

}

TEST(final, final_destination) {
	Sampler sampler(640, 480);
	Point cameraPos(0, 0, 500);
	Point lookingAt(0, 0, 300);
	Vector up(0, 1, 0);
	float fov = pi/2;
	Camera camera(640, 480, cameraPos, lookingAt, up, fov);
	TransMatrix mat;
	Point cen(-1, 0, 0);
	mat.add_translation(0, 390, -175);
	// mat.add_rotation(cen, pi/2);
	// mat.add_scaling(2, 1, 1);
	Transformation t(mat);
	Color kd(.8, .8, .3);
	Color ks(.9, .4, .6);
	Color ka(.4, .2, .4);
	Color kr(.2, .3, .3);
	float sp = 8;
	BRDF brdf(kd, ks, ka, kr, sp);
	Material* material = new Material(brdf);
	Color kd2(.7, .7, .3);
	Color ks2(.8, .6, .8);
	Color ka2(.1, .2, .2);
	Color kr2(.3, .3, .3);
	float sp2 = 16;
	BRDF brdf2(kd2, ks2, ka2, kr2, sp2);
	Material* material2 = new Material(brdf2);
	Point center1(0, -200.0, 0);
	Point center2(0, 0, 0);
	Point vertex1(-300, 200, -400);
	Point vertex2(-600, -300, -100);
	Point vertex3(0, -400, -300);
	Triangle* triangle = new Triangle(vertex1, vertex2, vertex3);
	GeometricPrimitive* gamma = new GeometricPrimitive(t, triangle, material);

	Point vertex4(300, 400, 400);
	Point vertex5(600, -200, 100);
	Point vertex6(0, -200, 200);
	Triangle* triangle2 = new Triangle(vertex4, vertex5, vertex6);
	Point center3(200, 200, 200);
	float radius3 = 50;
	float radius = 100;
	float radius2 = 150;
	Sphere* sphere1 = new Sphere(center1, radius);
	Sphere* sphere2 = new Sphere(center2, radius2);
	Sphere* sphere3 = new Sphere(center3, radius);
	GeometricPrimitive* alpha = new GeometricPrimitive(t, sphere3, material);
	GeometricPrimitive* beta = new GeometricPrimitive(t, sphere2, material2);
	GeometricPrimitive* delta = new GeometricPrimitive(t, triangle2, material2);
	vector<Primitive*> primitives;
	primitives.push_back(alpha);
	primitives.push_back(beta);
	primitives.push_back(gamma);
	primitives.push_back(delta);
	HBB prims(primitives, 0);
	Point light1 (500, 500, 500);
	Point light2 (0, -1000, 0);
	Color l1(1, 1, 1);
	Color l2(.6, .6, .6);
	PointLight* lig = new PointLight(light1, l1);
	PointLight* lig2 = new PointLight(light2, l2);
	vector<Light*> lights;
	lights.push_back(lig);
	lights.push_back(lig2);
	RayTracer raytracer(prims, lights, cameraPos);
	CImg<float> img(640, 480, 1, 3);
	const char* test = "tester2.jpg";
	Film film(img, false, test);
	Scene scene(sampler, camera, raytracer, film, 5);
	scene.render();
}

int main( int argc, char* argv[] )
{
	google::ParseCommandLineFlags(&argc, &argv, true);
	google::InitGoogleLogging(argv[0]);
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();

}
