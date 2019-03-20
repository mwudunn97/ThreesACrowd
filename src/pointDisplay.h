//
// Created by Marc WuDunn on 3/17/19.
//

#ifndef THREESACROWD_POINTDISPLAY_H
#define THREESACROWD_POINTDISPLAY_H

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <GL/freeglut.h>
#include <iterator>
#include <fstream>

#include <glm/vec3.hpp> // glm::vec3
#include <vector>
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale, glm::perspective
#include "Person.h"
using namespace glm;
using namespace std;


void display_point(glm::vec2 point);
int write_points(vector<vector<vec2>> points);
void draw_points(std::vector<std::vector<vec2>> points);
void set_points_from_groups(std::vector<Group> groups);
vector<vector<vec2>> generate_examples();
void display();
void myinit(float width, float height);
int pointDisplay(int argc, char** argv);
void display_points(float width, float height);
vector<vec2> points_from_groups(vector<Group> groups);
void set_points(vector<vector<vec2>> points);

#endif //THREESACROWD_POINTDISPLAY_H
