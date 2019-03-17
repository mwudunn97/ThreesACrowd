////
////  pointDisplay.cpp
////  ThreesACrowd
////
////  Created by Marc WuDunn on 2/23/19.
////  Copyright © 2019 Marc WuDunn. All rights reserved.
////


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <GLUT/glut.h>
#include <GLFW/glfw3.h>
#include <iterator>
#include <fstream>

#include <glm/vec3.hpp> // glm::vec3
#include <vector>
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale, glm::perspective
using namespace glm;
using namespace std;


const int num_circle_segments = 100;
const int circle_display_rad = 5.0;

void display_point(glm::vec2 point) {


        for(int ii = 0; ii < num_circle_segments; ii++)
        {
            float theta = 2.0f * 3.1415926f * float(ii) / float(num_circle_segments);//get the current angle

            float x = circle_display_rad * cosf(theta) + point[0];//calculate the x component
            float y = circle_display_rad * sinf(theta) + point[1];//calculate the y component

            glVertex2f(x, y);//output vertex

        }

}

int write_points(vector<vector<vec2>> points) {
    ofstream pointFile;
    pointFile.open ("locations.txt");
    pointFile << to_string(points[0].size());
    for (vector<vec2> trajectories: points) {
        for (vec2 point: trajectories) {
            pointFile << to_string(point[0]) + " " + to_string(point[1]) + "\n";
        }
    }
    pointFile.close();
    return 0;
}

void draw_points(std::vector<vec2> points) {
    GLfloat point_color[3] = {0.0,0.0,0.0};
    glEnable( GL_POINT_SMOOTH );
    glEnable( GL_BLEND );
    glPointSize(circle_display_rad);

    glColor3fv(point_color);
    glBegin(GL_POINTS);
        for (std::vector<vec2>::iterator vert = points.begin(); vert != points.end(); ++vert) {
        display_point(*vert);
    }
    glEnd();
    glFlush();
}

vector<vector<vec2>> generate_examples() {
    vector<vector<vec2>> traj;
    vector<vec2> ex_vec;
    ex_vec.push_back(glm::vec2(220.0, 222.0));
    ex_vec.push_back(glm::vec2(320.0, 200.0));
    ex_vec.push_back(glm::vec2(100.0, 325.0));
    ex_vec.push_back(glm::vec2(400.0, 300.0));
    traj.push_back(ex_vec);

    vector<vec2> ex_vec2;
    ex_vec2.push_back(glm::vec2(225.0, 222.0));
    ex_vec2.push_back(glm::vec2(325.0, 200.0));
    ex_vec2.push_back(glm::vec2(105.0, 325.0));
    ex_vec2.push_back(glm::vec2(405.0, 300.0));
    traj.push_back(ex_vec2);

    vector<vec2> ex_vec3;
    ex_vec3.push_back(glm::vec2(235.0, 222.0));
    ex_vec3.push_back(glm::vec2(335.0, 200.0));
    ex_vec3.push_back(glm::vec2(115.0, 325.0));
    ex_vec3.push_back(glm::vec2(415.0, 300.0));
    traj.push_back(ex_vec3);

    return traj;
}
void display()
{

    vector<vector<vec2>> example_traj = generate_examples();
    vector<vec2> vec = example_traj[0];

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBegin(GL_POINTS);
    draw_points(vec);

    glEnd();
    glFlush();

}

void myinit() {
    glEnable(GL_DEPTH_TEST);
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glColor3f(0.0, 1.0, 0.0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0,800,0.0,800.0);

    glMatrixMode(GL_MODELVIEW);
}

int pointDisplay(int argc, char** argv)
{

    vector<vector<vec2>> exam_points = generate_examples();
    write_points(exam_points);

    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);

    glutInitWindowSize(800, 800);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("GLUT Program");

    glutDisplayFunc(display);

    myinit();
    glutMainLoop();
    return 0;
}
