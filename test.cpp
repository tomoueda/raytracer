#include <iostream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include "CImg.h"
using namespace cimg_library;

DEFINE_bool(big_menu, true, "Include 'advanced' options in the menu listing");
DEFINE_string(languages, "english,french,german", "comma-seperated list of languages to offer in the 'lang' menu");

TEST(PrimitiveTest, primitive) {
	EXPECT_EQ(1, 1);
	EXPECT_EQ(1, 1);
}

int main( int argc, char* argv[] )
{
	google::ParseCommandLineFlags(&argc, &argv, true);
	google::InitGoogleLogging(argv[0]);
	LOG(INFO) << FLAGS_big_menu;
	// CImg<unsigned char> image("lena.jpg"), visu(500, 400, 1, 3, 0);
	// const unsigned char red[] = {255, 0, 0}, green[] = {0, 255, 0}, blue[] = {0, 0, 255};
	// image.blur(2.5);
	// CImgDisplay main_disp(image, "Click a point"), draw_disp(visu, "Intensity profile");
	// while (!main_disp.is_closed() && !draw_disp.is_closed()) {
	// 	main_disp.wait();
	// 	if(main_disp.button() && main_disp.mouse_y() >= 0) {
	// 		const int y = main_disp.mouse_y();
	// 		visu.fill(0).draw_graph(image.get_crop(0, y, 0, 0, image.width()-1, y, 0, 0), red, 1, 1, 0, 255, 0);
	// 		visu.draw_graph(image.get_crop(0, y, 0, 1, image.width()-1, y, 0, 1), green, 1, 1, 0, 255, 0);
	// 		visu.draw_graph(image.get_crop(0, y, 0, 2, image.width()-1, y, 0, 2), blue, 1, 1, 0, 255, 0).display(draw_disp);
	// 	}
	// }
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();

}