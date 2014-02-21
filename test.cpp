#include <iostream>
#include <gflags/gflags.h>

DEFINE_bool(big_menu, true, "Include 'advanced' options in the menu listing");
DEFINE_string(languages, "english,french,german", "comma-seperated list of languages to offer in the 'lang' menu");

int main( int argc, const char* argv[] )
{
	printf( "%d, \n", FLAGS_big_menu);
}