#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h> 
#include <stdio.h> 
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include <GraphicsMagick/Magick++.h>
#include <CImg.h>
#include "as2.h"

using namespace std;
/*
    Static Variable
*/
#define EPSILON .001

/*
    Flags
*/
DEFINE_string(input, "null", "an input file with the extension .txt");
DEFINE_bool(unittest, false, "set equal to one in order to run unit test.");

/* 
    Global Variables
*/
/** Image size width. **/
unsigned int width;
/** Image size height.**/
unsigned int height;
/** The maximum number of bounces for the ray **/
unsigned int depth;
/** Output file name. **/
std::string filename;
/** The camera origin. **/
float lookfromx, lookfromy, lookfromz;
/** The direction we are looking at. **/
float lookatx, lookaty, lookatz;
/** DON'T REALLY KNOW WHAT THIS IS FOR. **/
float upx, upy, upz;
/** Field of View. **/
float fov;

 

/*
    Utilitiy Functions.
*/
/** A function that sets up the function stream, and
    parses the inputs. **/
void parseInput() {
    if (FLAGS_input.compare("null")) {
        ifstream file(FLAGS_input.c_str());
        LOG(INFO) << "We are here";
        if (!file) {
            printf("%s does not exist, please make sure you're file name is correct", FLAGS_input.c_str());
            exit(EXIT_FAILURE);
        } else {
            LOG(INFO) << "We are here";
            string line;
            vector<string> parsedline;
            while (getline(file, line)) {
                if (line.compare("")) {
                    parsedline = split(line, ' ');
                    if (parsedline.at(0).compare("#")) {
                       parseInputLine(parsedline);
                    }
                }
            }
            LOG(INFO) << "WE NO HERE";
        }
    } else {
        printf("There were no file input, please input a file.\n");
        exit(EXIT_FAILURE);
    }
}

/** Takes in LINE which then parses according to the input. **/
void parseInputLine(std::vector<std::string> line) {
    LOG(INFO) << line.at(0);
    if (!line.at(0).compare("size")) {
        width = atoi(line.at(1).c_str());
        height = atoi(line.at(2).c_str());
    } else if (!line.at(0).compare("maxdepth")) {
        depth = atoi(line.at(1).c_str());
    } else if (!line.at(0).compare("output")) {
        filename = line.at(1);
        LOG(INFO) << "still works?";
    } else if (!line.at(0).compare("camera")) {
        LOG(INFO) << "so we work here?";
        lookfromx = atof(line.at(1).c_str());
        lookfromy = atof(line.at(2).c_str());
        lookfromz = atof(line.at(3).c_str());
        lookatx = atof(line.at(4).c_str());
        lookaty = atof(line.at(5).c_str());
        lookatz = atof(line.at(6).c_str());
        upx = atof(line.at(7).c_str());
        upy = atof(line.at(8).c_str());
        upz = atof(line.at(9).c_str());
        fov = atof(line.at(10).c_str());
    }
}

/** A function that splits STRING, by a DELIMiter, and stores it 
    in ELEM. Returns elem, which is a vector containing the split
    string. Disclaimer found this method on stack overflow. **/
std::vector<std::string> &split(const string &s, char delim, 
    std::vector<std::string> &elems) {
    stringstream sstream(s);
    string item;
    while (getline(sstream, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

/** Same as the following but a function that returns the a new Vector.
    **/
std::vector<std::string> split(const string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


int main(int argc, char **argv) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    if (FLAGS_unittest) {
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    }
    return 0;
}





/** UNIT TESTING. **/
TEST(utiltesting, split) {
    std::vector<std::string> test;
    std::string edge1 = " ";
    std::string edge2 = " ";
    test = split(edge1, ' ');
    LOG(INFO) << test.at(0) << "Are we here?";
    LOG(INFO) << " ";
    EXPECT_EQ(0, strcmp(" ", &edge1.at(0))) << "not empty string";
    test = split(edge2, ' ');
    LOG(INFO) << test.at(0) << "ARe WE HERE?";
    LOG(INFO) << " ";
    EXPECT_EQ(0, strcmp(" ", &edge2.at(0))) << "not empty string 2";
}

TEST(test1, initialParsingCheck) {
    FLAGS_input = "test1.txt";
    parseInput();
    LOG(INFO) << "We break here.";
    EXPECT_EQ(width, 400) << "Width is not the right value.";
    EXPECT_EQ(height, 400);
    EXPECT_EQ(depth, 5);
    EXPECT_EQ(lookfromx, 200);
    EXPECT_EQ(lookfromy, 100);
    EXPECT_EQ(lookfromz, 200);
    EXPECT_EQ(lookatx, 400);
    EXPECT_EQ(lookaty, 200);
    EXPECT_EQ(lookatz, 300);
    EXPECT_TRUE(fabs(upx - 1.1) < EPSILON);
    EXPECT_TRUE(fabs(upy - 2.1) < EPSILON);
    EXPECT_TRUE(fabs(upz - 3.3) < EPSILON);
    EXPECT_TRUE(fabs(fov - 18.9) < EPSILON) << "Is this roundoff error?";
}