//
// Created by Marc WuDunn on 3/17/19.
//

#ifndef THREESACROWD_POINTDISPLAY_H
#define THREESACROWD_POINTDISPLAY_H

#include <iostream>
#include <stdlib.h>
#include <math.h>
#ifdef __MINGW32__
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif
#include <iterator>
#include <fstream>

#include <glm/vec3.hpp> // glm::vec3
#include <vector>
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale, glm::perspective
using namespace glm;
using namespace std;


void display_point(glm::vec2 point);
int write_points(vector<vector<vec2>> points);
void draw_points(std::vector<vec2> points);
vector<vector<vec2>> generate_examples();
void display();
void myinit();
int pointDisplay(int argc, char** argv);
void display_points();

#endif //THREESACROWD_POINTDISPLAY_H
