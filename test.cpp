#include <iostream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include "CImg.h"
#include "classes.h"
#include "as2.h"

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
	a.debug();
	d.debug();
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
}

TEST(classes, normal) {
	Vector a(.4, 0, 0);
	Normal q(.1230921, 0, 0);
	Normal b(1, 0, 0);
	Normal c = a.normalize();
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
	a.debug();
	Vector y = a - b;
	a.debug();
	EXPECT_EQ(y, z) << "Point - Point = Vector";
	a = a + z;
	Point c(-8, -6, -4);
	a.debug();
	c.debug();
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
	a.debug();
	d.debug();
	EXPECT_EQ(d, c) << "testing '+' operator";
	EXPECT_EQ(a, e) << "Making sure that a is untouched. ";
	a += b;
	a.debug();
	d.debug();
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
	BRDF brdf(kd, ks, ka, kr);
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

int main( int argc, char* argv[] )
{
	google::ParseCommandLineFlags(&argc, &argv, true);
	google::InitGoogleLogging(argv[0]);
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();

}