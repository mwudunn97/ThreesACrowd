//
// Created by Dalton on 2/24/2019.
//

#ifndef THREESACROWD_GRID_H
#define THREESACROWD_GRID_H

#include <glm/vec4.hpp>
#include <glm/vec2.hpp>
#include <vector>
#include <array>
#include <unordered_map>
#include <json.hpp>
#include "Person.h"

using json = nlohmann::json;

enum Direction {
  East = 0,
  North = 1,
  West = 2,
  South = 3
};

struct Edge {
  /* Positive is Northward and Eastward */
  float h_grad;
  float phi_grad;
  float v;
};

/* n sub theta vector representing unit directions in the
 * order in the paper: E, N, W, S.
 * Direction "OUT" of a cell */
const std::array<int, 4> n_theta_int {1, 1, -1, -1};
const std::array<glm::vec2, 4> n_theta_vec {{
  {1, 0}, {0, 1}, {-1, 0}, {0, -1}
}};

struct Cell {
  Cell(Edge *edgeE, Edge *edgeN, Edge *edgeW, Edge *edgeS, int i, int j);

  /* Location (for debugging) */
  int i;
  int j;

  /* Center of this cell */
  float g = 0;
  float phi = 0;
  float rho = 0;
  float h = 0;
  glm::vec2 v_avg {};

  /* Going OUT OF this cell: E,N,W,S */
  glm::vec4 f {};
  glm::vec4 C {};

  std::array<Edge*, 4> edges;

  /* Neighbors */
  std::array<Cell*, 4> neighbors;
};

/* Row-major grid
 * 0,0 in bottom left */
class Grid {
public:
  Grid(int width, int height);
  explicit Grid(json &j);

  int getWidth() const;
  int getHeight() const;

  void clearGridVals();

  Cell *getCell(int i, int j);
  Cell *getCell(glm::ivec2 ij);

  std::vector<std::vector<Cell>> grid;
  std::vector<Edge> edges;
  std::unordered_map<float, std::vector<Person *> *> map;

  int width;
  int height;

  /* Tunable variables */
  float alpha;   // eqn. 4
  float beta;    // eqn. 4
  float gamma;   // eqn. 4
  float rho_min; // eqn. 10
  float rho_max; // eqn. 10
  float f_min;   // eqn. 8
  float f_max;   // eqn. 8
  float s_min;   // eqn. 8
  float s_max;   // eqn. 8
  double lambda; // section 4.1

  void build_neighbor_map(std::vector<Group> &groups);
  void handle_collisions(Person &person);

private:
  void fill();
  float hash_position(glm::vec2 pos);
};


#endif //THREESACROWD_GRID_H
