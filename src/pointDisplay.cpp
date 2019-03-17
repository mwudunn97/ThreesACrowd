////
////  pointDisplay.cpp
////  ThreesACrowd
////
////  Created by Marc WuDunn on 2/23/19.
////  Copyright © 2019 Marc WuDunn. All rights reserved.
////


#include "pointDisplay.h"


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

//Write points to a file, so it can be loaded into unity
int write_points(vector<vector<vec2>> points, string filename) {
  ofstream pointFile;
  pointFile.open (filename);
  pointFile << to_string(points[0].size());
  for (vector<vec2> trajectories: points) {
    for (vec2 point: trajectories) {
      pointFile << to_string(point[0]) + " " + to_string(point[1]) + "\n";
    }
  }
  pointFile.close();
  return 0;
}

//Takes a set of input points and draws them to the screen
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

//Test example
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

//GLUT needs a global vec
vector<vector<vec2>> point_traj;

void set_points(vector<vector<vec2>> points) {
  point_traj = points;
}

vector<vector<vec2>> get_points() {
  return point_traj;
}

//GLUT display function
void display()
{

  vector<vector<vec2>> points = get_points();
  vector<vec2> vec = points[0];

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glBegin(GL_POINTS);
  draw_points(vec);

  glEnd();
  glFlush();

}

//Init matrices
void myinit() {
  glEnable(GL_DEPTH_TEST);
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glColor3f(0.0, 1.0, 0.0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0,800,0.0,800.0);

  glMatrixMode(GL_MODELVIEW);
}

void display_points(vector<vector<vec2>> points) {
  set_points(points);
  glutDisplayFunc(display);

  myinit();
  glutMainLoop();

}

//Call this function to display points
int pointDisplay(int argc, char** argv)
{

  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);

  glutInitWindowSize(800, 800);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("GLUT Program");
  return 0;

}
