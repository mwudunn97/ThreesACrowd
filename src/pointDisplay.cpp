////
////  pointDisplay.cpp
////  ThreesACrowd
////
////  Created by Marc WuDunn on 2/23/19.
////  Copyright © 2019 Marc WuDunn. All rights reserved.
////


#include "pointDisplay.h"


const int num_circle_segments = 1000;
const int circle_display_rad = 1.0;
const float resX = 100.0f;
const float resY = 100.0f;

//GLUT needs a global vec
vector<vector<vec2>> point_traj;
float gridWidth = 10.0f;
float gridHeight = 10.0f;

void display_point(glm::vec2 point) {


  for(int ii = 0; ii < num_circle_segments; ii++)
  {
    float theta = 2.0f * 3.1415926f * float(ii) / float(num_circle_segments);//get the current angle

    float x = circle_display_rad * cosf(theta) + point[0];//calculate the x component
    float y = circle_display_rad * sinf(theta) + point[1];//calculate the y component

    glVertex2f(x, y);//output vertex

  }

}

void set_grid_sizes(float width, float height) {
  gridWidth = width;
  gridHeight = height;
}

//Write points to a file, so it can be loaded into unity
int write_points(vector<vector<vec2>> points, string filename) {
  ofstream pointFile;
  pointFile.open (filename);
  std::cout << points[0].size();
  pointFile << to_string(points[0].size()) + " " + to_string(points.size()) + "\n";
  for (vector<vec2> trajectories: points) {
    for (vec2 point: trajectories) {
      pointFile << to_string(point[0] / gridWidth) + " " + to_string(point[1] / gridHeight) + "\n";
    }
  }
  pointFile.close();
  return 0;
}

//Takes a set of input points and draws them to the screen
void draw_points(vector<vector<vec2>> point_traj) {
  GLfloat point_color[3] = {0.0,0.0,0.0};
  glEnable( GL_POINT_SMOOTH );
  glEnable( GL_BLEND );
  glPointSize(circle_display_rad);

  glColor3fv(point_color);
  glBegin(GL_POINTS);
  for (auto &points : point_traj) {
    for (std::vector<vec2>::iterator vert = points.begin(); vert != points.end(); ++vert) {
      vec2 displayVert = glm::vec2(vert->x * (resX / gridWidth), vert->y * (resY / gridHeight));
      display_point(displayVert);
    }
  }
  glEnd();
  glFlush();
}

/* Draws people as points */
void draw_people(std::vector<Group> groups) {
  GLfloat point_color[3] = {0.0,0.0,0.0};
  glEnable( GL_POINT_SMOOTH );
  glEnable( GL_BLEND );
  glPointSize(circle_display_rad);

  glColor3fv(point_color);
  glBegin(GL_POINTS);
  for (auto &group : groups) {
    for (auto &person : group.people) {
      display_point(person.getPos());
    }
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





//Use this to set the point vector that gets drawn
void set_points(vector<vector<vec2>> points) {
  point_traj = points;
}

//Alternative function that operates on groups
void set_points_from_groups(vector<Group> groups) {

  vector<vec2> people_pos;
  for (auto &group : groups) {
    for (auto &person : group.people) {
      people_pos.push_back(person.getPos());
    }
  }

  vector<vector<vec2>> points;
  points.push_back(people_pos);
  point_traj = points;

}

//Converts groups to single vector of points
vector<vec2> points_from_groups(vector<Group> groups) {

  vector<vec2> people_pos;
  for (auto &group : groups) {
    for (auto &person : group.people) {
      people_pos.push_back(person.getPos());
    }
  }

  return people_pos;

}



vector<vector<vec2>> get_points() {
  return point_traj;
}

//GLUT display function, anything you want to display has to be called in this function
void display()
{

  vector<vector<vec2>> points = get_points();

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glBegin(GL_POINTS);
  draw_points(points);

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
  gluOrtho2D(0.0,resX,0.0,resY);

  glMatrixMode(GL_MODELVIEW);
}

void display_points(float width, float height) {
  glutDisplayFunc(display);
  gridWidth = width;
  gridHeight = height;
  myinit();
  /* Allows you to call the display function and continue processing, only works if freeglut works on your machine */
  //glutMainLoopEvent();
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
