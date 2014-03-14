
#include <time.h>
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
#include "classes.h"

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
unsigned int width = 640;
/** Image size height.**/
unsigned int height = 480;
/** Our sampler to generate sample points. **/
Sampler sampler;
/** The maximum number of bounces for the ray **/
unsigned int depth = 5;
/** Output file name. **/
string filename;
/** The camera used for the scene. **/
Camera camera;
/** A list of primitives. **/
vector<Primitive*> primitives;
/** A list of vertex, stored as points. **/
vector<Point> vertexes;
/** A list of normalVertexes, stored as vertex. **/
vector<Vertex> Nvertexes;
/** A list of normals. **/
vector<Normal> normals;
/** A boolean to check whether the next transformation will be added. **/
bool canPush = false;
/** Our current TransMatrix **/
TransMatrix currMatrix;
/** Lights stored in a vector. **/
vector<Light*> lights;
/** The current colors. **/
Color kd, ks, ka, kr;
/** The shininess. **/
float sp;

 

/*
    Utilitiy Functions.
*/
/** A function that sets up the function stream, and
    parses the inputs. **/
void parseInput() {
    if (FLAGS_input.compare("null")) {
        ifstream file(FLAGS_input.c_str());
        if (!file) {
            printf("%s does not exist, please make sure you're file name is correct", FLAGS_input.c_str());
            exit(EXIT_FAILURE);
        } else {
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
        }
    } else {
        printf("There were no file input, please input a file.\n");
        exit(EXIT_FAILURE);
    }
}

/** Takes in LINE which then parses according to the input. **/
void parseInputLine(std::vector<std::string> line) {
    if (!line.at(0).compare("size")) {
        width = atoi(line.at(1).c_str());
        height = atoi(line.at(2).c_str());
        sampler = *new Sampler(width, height);
    } else if (!line.at(0).compare("maxdepth")) {
        depth = atoi(line.at(1).c_str());
    } else if (!line.at(0).compare("output")) {
        filename = line.at(1);
        LOG(INFO) << filename;
    } else if (!line.at(0).compare("camera")) {
        float lookfromx = atof(line.at(1).c_str());
        float lookfromy = atof(line.at(2).c_str());
        float lookfromz = atof(line.at(3).c_str());
        Point lookfrom(lookfromx, lookfromy, lookfromz);
        float lookatx = atof(line.at(4).c_str());
        float lookaty = atof(line.at(5).c_str());
        float lookatz = atof(line.at(6).c_str());
        Point lookat(lookatx, lookaty, lookatz);
        float upx = atof(line.at(7).c_str());
        float upy = atof(line.at(8).c_str());
        float upz = atof(line.at(9).c_str());
        Vector up(upx, upy, upz);
        float fov = atof(line.at(10).c_str());
        camera = *new Camera(width, height, lookfrom, lookat, up, fov);
    } else if (!line.at(0).compare("sphere")) {
        float x = atof(line.at(1).c_str());
        float y = atof(line.at(2).c_str());
        float z = atof(line.at(3).c_str());
        float radius = atof(line.at(4).c_str());
        Point center(x, y, z);
        Sphere* sphere = new Sphere(center, radius);
        Transformation t(currMatrix);
        BRDF brdf(kd, ks, ka, kr, sp);
        Material* material = new Material(brdf);
        GeometricPrimitive* s = new GeometricPrimitive(t, sphere, material);
        primitives.push_back(s);
    } else if (!line.at(0).compare("vertex") || !line.at(0).compare("v")) {
        float x = atof(line.at(1).c_str());
        float y = atof(line.at(2).c_str());
        float z = atof(line.at(3).c_str());
        Point vertex(x, y, z);
        vertexes.push_back(vertex);
    } else if (!line.at(0).compare("vn")) {
        float x = atof(line.at(1).c_str());
        float y = atof(line.at(2).c_str());
        float z = atof(line.at(3).c_str());
        Normal normal(x, y, z);
        normals.push_back(normal);
    } else if (!line.at(0).compare("f")) {
        if (line.at(1).find('/') != string::npos) {
            vector<string> sv1 = split(line.at(1), '/');
            vector<string> sv2 = split(line.at(2), '/');
            vector<string> sv3 = split(line.at(3), '/');
            Vertex v1(vertexes.at(atoi(sv1.at(0).c_str()) - 1), normals.at(atoi(sv1.at(2).c_str()) - 1));
            Vertex v2(vertexes.at(atoi(sv2.at(0).c_str()) - 1), normals.at(atoi(sv2.at(2).c_str()) - 1));
            Vertex v3(vertexes.at(atoi(sv3.at(0).c_str()) - 1), normals.at(atoi(sv3.at(2).c_str()) - 1));
            VertexTriangle* tri = new VertexTriangle(v1, v2, v3);
            Transformation t(currMatrix);
            BRDF brdf(kd, ks, ka, kr, sp);
            Material* material = new Material(brdf);
            GeometricPrimitive* vertri = new GeometricPrimitive(t, tri, material);
            primitives.push_back(vertri);
        } else {
            int len = line.size();
            if (len == 4) {
                Point v1 = vertexes.at(atoi(line.at(1).c_str())-1);
                Point v2 = vertexes.at(atoi(line.at(2).c_str())-1);
                Point v3 = vertexes.at(atoi(line.at(3).c_str())-1);
                Triangle* tri = new Triangle(v1, v2, v3);
                Transformation t(currMatrix);
                BRDF brdf(kd, ks, ka, kr, sp);
                Material* material = new Material(brdf);
                GeometricPrimitive* triangle = new GeometricPrimitive(t, tri, material);
                primitives.push_back(triangle);
            } else if (len == 5) {
                Point v1 = vertexes.at(atoi(line.at(1).c_str())-1);
                Point v2 = vertexes.at(atoi(line.at(2).c_str())-1);
                Point v3 = vertexes.at(atoi(line.at(3).c_str())-1);
                Triangle* tri = new Triangle(v1, v2, v3);
                Transformation t(currMatrix);
                BRDF brdf(kd, ks, ka, kr, sp);
                Material* material = new Material(brdf);
                GeometricPrimitive* triangle = new GeometricPrimitive(t, tri, material);
                primitives.push_back(triangle);

                Point v4 = vertexes.at(atoi(line.at(2).c_str())-1);
                Point v5 = vertexes.at(atoi(line.at(3).c_str())-1);
                Point v6 = vertexes.at(atoi(line.at(4).c_str())-1);
                Triangle* tri2 = new Triangle(v4, v5, v6);
                Transformation t2(currMatrix);
                BRDF brdf2(kd, ks, ka, kr, sp);
                Material* material2 = new Material(brdf2);
                GeometricPrimitive* triangle2 = new GeometricPrimitive(t2, tri2, material2);
                primitives.push_back(triangle2);
            }
        }
    } else if (!line.at(0).compare("vertexnormal")) {
        float x = atof(line.at(1).c_str());
        float y = atof(line.at(2).c_str());
        float z = atof(line.at(3).c_str());
        float nx = atof(line.at(4).c_str());
        float ny = atof(line.at(5).c_str());
        float nz = atof(line.at(6).c_str());
        Point pos(x, y, z);
        Normal norm(nx, ny, nz);
        Vertex vex(pos, norm);
        Nvertexes.push_back(vex);
    } else if (!line.at(0).compare("trinormal")) {
        Vertex v1 = Nvertexes.at(atoi(line.at(1).c_str()) - 1);
        Vertex v2 = Nvertexes.at(atoi(line.at(2).c_str()) - 1);
        Vertex v3 = Nvertexes.at(atoi(line.at(3).c_str()) - 1);
        VertexTriangle* vt = new VertexTriangle(v1, v2, v3);
        Transformation t(currMatrix);
        BRDF brdf(kd, ks, ka, kr, sp);
        Material* material = new Material(brdf);
        GeometricPrimitive* vtgp = new GeometricPrimitive(t, vt, material);
        primitives.push_back(vtgp);
    } else if (!line.at(0).compare("tri")) {
        Point v1 = vertexes.at(atoi(line.at(1).c_str())-1);
        Point v2 = vertexes.at(atoi(line.at(2).c_str())-1);
        Point v3 = vertexes.at(atoi(line.at(3).c_str())-1);
        Triangle* tri = new Triangle(v1, v2, v3);
        Transformation t(currMatrix);
        BRDF brdf(kd, ks, ka, kr, sp);
        Material* material = new Material(brdf);
        GeometricPrimitive* triangle = new GeometricPrimitive(t, tri, material);
        primitives.push_back(triangle);
    } else if (!line.at(0).compare("pushTransform")) {
        canPush = true;
    } else if (!line.at(0).compare("popTransform")) {
        currMatrix = *new TransMatrix();
        canPush = false;
    } else if (!line.at(0).compare("translate")) {
        if (canPush) {
            float x = atof(line.at(1).c_str());
            float y = atof(line.at(2).c_str());
            float z = atof(line.at(3).c_str());
            currMatrix.add_translation(x, y, z);
        }
    } else if (!line.at(0).compare("rotate")) {
        if (canPush) {
            float x = atof(line.at(1).c_str());
            float y = atof(line.at(2).c_str());
            float z = atof(line.at(3).c_str());
            float theta = atof(line.at(4).c_str());
            Point about(x, y, z);
            currMatrix.add_rotation(about, theta);
        }
    } else if (!line.at(0).compare("scale")) {
        if (canPush) {
            float x = atof(line.at(1).c_str());
            float y = atof(line.at(2).c_str());
            float z = atof(line.at(3).c_str());
            currMatrix.add_scaling(x, y, z);
        }
    } else if (!line.at(0).compare("directional")) {
        float x = atof(line.at(1).c_str());
        float y = atof(line.at(2).c_str());
        float z = atof(line.at(3).c_str());
        float r = atof(line.at(4).c_str());
        float g = atof(line.at(5).c_str());
        float b = atof(line.at(6).c_str());
        Vector dir(x, y, z);
        Color intense(r, g, b);
        DirectionalLight* light = new DirectionalLight(dir, intense);
        lights.push_back(light);
    } else if (!line.at(0).compare("point")) {
        float x = atof(line.at(1).c_str());
        float y = atof(line.at(2).c_str());
        float z = atof(line.at(3).c_str());
        float r = atof(line.at(4).c_str());
        float g = atof(line.at(5).c_str());
        float b = atof(line.at(6).c_str());
        Point pos(x, y, z);
        Color intense(r, g, b);
        PointLight* light = new PointLight(pos, intense);
        lights.push_back(light);
    } else if (!line.at(0).compare("diffuse")) {
        float r = atof(line.at(1).c_str());
        float g = atof(line.at(2).c_str());
        float b = atof(line.at(3).c_str());
        Color intense(r, g, b);
        kd.setValue(r, g, b);
    } else if (!line.at(0).compare("specular")) {
        float r = atof(line.at(1).c_str());
        float g = atof(line.at(2).c_str());
        float b = atof(line.at(3).c_str());
        Color intense(r, g, b);
        ks.setValue(r, g, b);
    } else if (!line.at(0).compare("reflect")) {
        float r = atof(line.at(1).c_str());
        float g = atof(line.at(2).c_str());
        float b = atof(line.at(3).c_str());
        Color intense(r, g, b);
        kr.setValue(r, g, b);
    } else if (!line.at(0).compare("ambient")) {
        float r = atof(line.at(1).c_str());
        float g = atof(line.at(2).c_str());
        float b = atof(line.at(3).c_str());
        Color intense(r, g, b);
        ka.setValue(r, g, b);
    } else if (!line.at(0).compare("shininess")) {
        sp = atof(line.at(1).c_str());
    }
}

/** Sets up all necessary components and starts rendering. **/
void render() {
    HBB agg(primitives, 0);
    RayTracer rayTracer(agg, lights, camera.getCameraPos());
    CImg<float> img(width, height, 1, 3);
    Film film(img, false, filename.c_str());
    Scene scene(sampler, camera, rayTracer, film, depth);
    scene.render();
}

void renderUsingList() {
    AggregatePrimitive agg(primitives);
    RayTracer rayTracer(agg, lights, camera.getCameraPos());
    CImg<float> img(width, height, 1, 3);
    Film film(img, false, filename.c_str());
    Scene scene(sampler, camera, rayTracer, film, depth);
    scene.renderUsingList();
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
    clock_t start, end;

    start = clock();
    parseInput();
    render();
    end = clock();
    LOG(INFO) << "HBB: Time required for execution: " << (double)(end-start)/CLOCKS_PER_SEC << " seconds.";

    // start = clock();
    // parseInput();
    // renderUsingList();
    // end = clock();
    // LOG(INFO) << "LIST: Time required for execution: " << (double)(end-start)/CLOCKS_PER_SEC << " seconds.";
    return 0;
}