/*
 * ThreesACrowd - implementation of "Continuum Crowds", Treuille, et al.
 * https://grail.cs.washington.edu/projects/crowd-flows/78-treuille.pdf
 * Dalton Omens, Marc WuDunn, Jessie Yang
*/

#include <glm/glm.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include "Grid.h"
#include "Person.h"
#include <json.hpp>
#include <GL/glut.h>
#include <algorithm>

using json = nlohmann::json;

void density_conversion(Grid &grid, Group &group, double lambda) {
  // TODO 4.1: convert positions of Persons into densities and
  // insert into Grid. Also calculate average velocities of each cell.

  std::vector<Person> &people = group.people;

  int width = grid.getWidth();
  int height = grid.getHeight();

  // reset all rho and v_avg values for the grid
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      Cell *cell = grid.getCell(i, j);
      cell->rho = 0.0f;
      cell->v_avg[0] = 0.0f;
      cell->v_avg[1] = 0.0f;
    }
  }


  for (auto &person : people) {
    // TODO: will this find the correct, closest cell center?
    glm::ivec2 gridIndex = person.getGridIndex();
    float dx = person.getPos()[0] - gridIndex[0];
    float dy = person.getPos()[1] - gridIndex[1];
    Cell *curr_cell = grid.getCell(gridIndex[0], gridIndex[1]);

    // add density to current cell
    float rho_a = static_cast<float>(pow(std::min(1-dx, 1-dy), lambda));
    curr_cell->rho += rho_a;
    // accumulate weighted density for avg velocity calculation
    curr_cell->v_avg += person.getVelocity() * rho_a;

    // TODO: think of a cleaner way to write this
    if (gridIndex[0] + 1 < width) {
      // add density to cell to the right
      float rho_b = static_cast<float>(pow(std::min(dx, 1-dy), lambda));
      curr_cell->neighbors[East]->rho += rho_b;

      if (gridIndex[1] - 1 >= 0) {
        // add density to cell above and to the right
        float rho_c = static_cast<float>(pow(std::min(dx, dy), lambda));
        curr_cell->neighbors[North]->neighbors[East]->rho += rho_c;
      }
    }
    if (gridIndex[1] - 1 >= 0) {
      // add density to cell above
      float rho_d = static_cast<float>(pow(std::min(1-dx, dy), lambda));
        curr_cell->neighbors[North]->rho += rho_d;
    }
  }

  // calculate the average velocity
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      Cell *cell = grid.getCell(i, j);
      cell->v_avg /= cell->rho;
    }
  }
}



void calculate_unit_cost(Grid &grid) {
  // TODO 4.2: iterate over each of 4 directions of each cell
  // and calculate f_{M->i} and {C_M->i} where i is in {N, E, S, W}

  /* To make things easier, we iterate through speed/cost going INTO a cell because you have
  * to use the fields from the cell you're going into, not the one you came from */
  for (auto &row : grid.grid) {
    for (auto &cell : row) {
      for (int dir = East; dir <= South; dir++) {
        /* dir = direction we're coming FROM. Therefore, n_theta is negated */
        Edge *e = cell.edges[dir];

        /* Equation 8: topographical speed */
        float slope = (glm::dot(e->h_grad, -n_theta[dir]) - grid.s_min) /
            (grid.s_max - grid.s_min);
        float topo_speed = grid.f_max + slope * (grid.f_min - grid.f_max);

        /* Equation 9: flow speed */
        float flow_speed = glm::dot(cell.v_avg, -n_theta[dir]);
        if (flow_speed < 0.0f) flow_speed = 0.0f;

        if (cell.rho <= grid.rho_min) {
          cell.f[dir] = topo_speed;
        } else if (cell.rho >= grid.rho_max) {
          cell.f[dir] = flow_speed;
        } else {
          /* Equation 10: linearly interpolate */
          float speed = topo_speed + (cell.rho - grid.rho_min) /
              (grid.rho_max - grid.rho_min) * (flow_speed - topo_speed);
          cell.f[dir] = speed;
        }

        /* Equation 4: Cost field */
        float cost = (grid.alpha * cell.f[dir] + grid.beta + grid.gamma * cell.g) / cell.f[dir];
        cell.C[dir] = cost;
      }
    }
  }
}

std::vector<Cell*> flatten(Grid &grid) {
  std::vector<Cell*> result;
  for (int j = 0; j < grid.width; j++) {
    for (int i = 0; i < grid.height; i++) {
      result.push_back(grid.getCell(i, j));
    }
  }
  return result;
}


void finite_differences_approx(Cell &cell) {
  float phi_mx;
  float phi_my;
  Direction d_mx;
  Direction d_my;

  //Check boundary cases
  if (cell.neighbors[West] == NULL) {
    phi_mx = cell.neighbors[East]->phi + cell.C[East];
    d_mx = East;
  } else if (cell.neighbors[East] == NULL) {
    phi_mx = cell.neighbors[East]->phi + cell.C[East];
    d_mx = West;
  } else {
    //Otherwise, choose minimum phi + index between west/east directions
    float phi_wx = cell.neighbors[West]->phi + cell.C[West];
    float phi_ex = cell.neighbors[East]->phi + cell.C[East];

    if (phi_wx < phi_ex) {
      phi_mx = phi_wx;
      d_mx = West;
    } else {
      phi_mx = phi_ex;
      d_mx = East;
    }
  }
  //Check boundary cases
  if (cell.neighbors[South] == NULL) {
    phi_my = cell.neighbors[North]->phi + cell.C[North];
    d_my = North;
  } else if (cell.neighbors[North] == NULL) {
    phi_my = cell.neighbors[South]->phi + cell.C[South];
    d_my = South;
  } else {
    //Otherwise, choose minimum phi + index between north/south directions
    float phi_ny = cell.neighbors[North]->phi + cell.C[North];
    float phi_sy = cell.neighbors[South]->phi + cell.C[South];

    if (phi_sy < phi_ny) {
      phi_my = phi_sy;
      d_my = South;
    } else {
      phi_my = phi_ny;
      d_my = North;
    }
  }

  //Set the different terms of the quadratic equation
  float c_mx = cell.C[d_mx];
  float c_my = cell.C[d_my];
  float phi_m;

  //If one of the terms is undefined, remove it from the quadratic equation
  //Otherwise, calculate the phi_m for the cell, and update velocities using speed field
  if (std::isinf(phi_mx)) {
    float det = c_mx;
    phi_m = phi_my + std::sqrt(det);
    cell.edges[d_my]->phi_grad = phi_m - phi_my;
    cell.edges[d_mx]->v = cell.edges[d_my]->phi_grad * (float) cell.f[d_my];
  } else if (std::isinf(phi_my)) {
    float det = c_my;
    phi_m = phi_mx + std::sqrt(det);
    cell.edges[d_mx]->phi_grad = phi_m - phi_mx;
    cell.edges[d_mx]->v = cell.edges[d_mx]->phi_grad * (float) cell.f[d_mx];
  } else {
    float a = c_mx + c_my;
    float b = 2.0f * (c_my * phi_mx + c_mx * phi_my);
    float c = (c_my * phi_mx * phi_mx) + (c_mx * phi_my * phi_my) - (c_mx * c_my);
    float det = b * b - 4.0f * a * c;
    phi_m = (-1.0f * b + std::sqrt(det)) / (2.0f * a);

    cell.edges[d_mx]->phi_grad = phi_m - phi_mx;
    cell.edges[d_mx]->v = cell.edges[d_mx]->phi_grad * (float) cell.f[d_mx];
    cell.edges[d_my]->phi_grad = phi_m - phi_my;
    cell.edges[d_mx]->v = cell.edges[d_my]->phi_grad * (float) cell.f[d_my];
  }
    cell.phi = phi_m;




}

//Compare function for the heap
bool cmp(const Cell * a, const Cell * b) {
  //overloaded < compare, see Cell
  //tbh could have just compared it directly here but wrote this after rip
  return *a < *b;
}


void construct_dynamic_potential_field(Grid &grid, glm::ivec2 goal) {
  // TODO 4.3: use fast-marching algorithm to calculate the potentials
  // and potential gradients (phi and del phi)
  // also calculate the velocity field of velocities v
  std::vector<Cell*> flattened_grid = flatten(grid);
  int goal_index = goal[1] * grid.getWidth() + goal[0];
  flattened_grid[goal_index]->phi = 0;
  std::make_heap(flattened_grid.begin(), flattened_grid.end(), cmp);
  for (int i = 0; i < flattened_grid.size(); i++) {
    Cell* next = flattened_grid.front();
    std::pop_heap (flattened_grid.begin(),flattened_grid.end()); flattened_grid.pop_back();

    for (Cell *c : next->neighbors) {
      if (c != nullptr) finite_differences_approx(*c);
    }
    finite_differences_approx(*next);
    std::make_heap(flattened_grid.begin(), flattened_grid.end() - i);
  }

  return;
}


void crowd_advection(Grid &grid, Group &group) {
  // TODO 4.4: update each person's position by interpolating into the vector
  // field
  std::vector<Person> &people = group.people;

  for (auto &person : people) {
    person.setPos(person.getPos() + person.getCell(grid)->v_avg);
  }
}


void enforce_minimum_distance(Grid &grid, Group &group) {
  // TODO 4.5: iterate over all pairs in a threshold distance and push people
  // apart symmetrically until min. distance is reached. may instead use a
  // neighbor grid instead of the vector of Persons.

  // generate spatial map
  std::vector<Person> &people = group.people;
  grid.build_neighbor_map(people);

  // run pair-wise minimum distance enforcement for each person
  for (auto &person : people) {
    grid.handle_collisions(person);
  }
}

void test_structures() {
  std::cout << "Hello World" << std::endl;

  Grid grid(4, 3);

  Person dalton(2.6f, 1.3f, 0, 0, -5.0f);
  Cell *daltonCell = dalton.getCell(grid); // 2,1
  daltonCell->g = 1234;
  daltonCell->edges[North]->v = glm::vec2(0.69, 0.420);
  Cell *daltonAbove = grid.getCell(dalton.getGridIndex() + glm::ivec2(0, 1));
  std::cout << daltonAbove->edges[South]->v[0] << " " << daltonAbove->edges[South]->v[1] << std::endl;
  std::cout << daltonAbove->neighbors[South]->g << std::endl;

  float smth = glm::dot(glm::vec2(7, 9), n_theta[South]);
  std::cout << "-9: " << smth << std::endl;
}

void test_potential_field() {
  Grid grid(4, 3);

  Person dalton(2.6f, 1.3f, 0.0f, 0.0f, -5.0f);
  Cell *daltonCell = dalton.getCell(grid); // 2,1
  daltonCell->g = 1234;
  daltonCell->edges[North]->v = glm::vec2(0.69, 0.420);
  Cell *daltonAbove = grid.getCell(dalton.getGridIndex() + glm::ivec2(0, 1));
  std::cout << daltonAbove->edges[South]->v[0] << " " << daltonAbove->edges[South]->v[1] << std::endl;
  std::cout << daltonAbove->neighbors[South]->g << std::endl;

  float smth = glm::dot(glm::vec2(7, 9), n_theta[South]);
  construct_dynamic_potential_field(grid, glm::ivec2(0,0));
}

int load_config(json &j, char *config) {
  std::ifstream f(config);
  if (!f.good()) {
    std::cerr << "Config file " << config << " does not exist" << std::endl;
    return -1;
  }
  f >> j;
  return 0;
}

int load_groups(json &j, std::vector<Group> *groups) {
  for (auto &group : j["groups"]) {
    std::vector<Person> people;
    for (auto &person : group["people"]) {
      people.emplace_back(person[0], person[1], person[2], person[3], person[4]);
    }
    glm::vec2 goal(group["goal"][0], group["goal"][1]);
    groups->emplace_back(goal, people);
  }

  return 0;
}

int main(int argc, char* argv[]) {
  std::cout << "Three's A Crowd Simulator" << std::endl;
  // test_structures();
  test_potential_field();

  if (argc < 2) {
    std::cerr << "Please specify a configuration file" << std::endl;
    return -1;
  }
  json j;
//  if (load_config(j, argv[1])) {
//    return -1;
//  }

  Grid grid(j);
  std::vector<Group> groups;
  if (load_groups(j, &groups)) {
    return -1;
  }

  int iterations = j["iterations"];
  for (int i = 0; i < iterations; i++) {

  }
  return 0;
}
