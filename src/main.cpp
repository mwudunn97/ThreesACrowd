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
#include <algorithm>
#include "pointDisplay.h"
#include <iomanip>

using json = nlohmann::json;

void density_conversion(Grid &grid, std::vector<Group> &groups) {
  /* convert positions of Persons into densities and
     insert into Grid. Also calculate average velocities of each cell. */

  double lambda = grid.lambda;
  for (auto &group : groups) {
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
      /* Simplifying: we just take the cell index of the one we're in, not the center */
      glm::ivec2 gridIndex = glm::floor(person.getPos());
      float dx = person.getPos()[0] - gridIndex[0];
      float dy = person.getPos()[1] - gridIndex[1];
      Cell *curr_cell = grid.getCell(gridIndex);

      /* add density to current cell */
      auto rho_a = static_cast<float>(pow(std::min(1 - dx, 1 - dy), lambda));
      curr_cell->rho += rho_a;
      // accumulate weighted density for avg velocity calculation
      curr_cell->v_avg += person.getVelocity() * rho_a;

      // TODO: think of a cleaner way to write this
      if (gridIndex[0] + 1 < width) {
        /* add density to cell to the right */
        auto rho_b = static_cast<float>(pow(std::min(dx, 1 - dy), lambda));
        curr_cell->neighbors[East]->rho += rho_b;

        if (gridIndex[1] + 1 < height) {
          /* add density to cell above and to the right */
          auto rho_c = static_cast<float>(pow(std::min(dx, dy), lambda));
          curr_cell->neighbors[North]->neighbors[East]->rho += rho_c;
        }
      }

      if (gridIndex[1] + 1 < height) {
        /* add density to cell above */
        auto rho_d = static_cast<float>(pow(std::min(1 - dx, dy), lambda));
        curr_cell->neighbors[North]->rho += rho_d;
      }
    }


    // calculate the average velocity
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        Cell *cell = grid.getCell(i, j);
        if (cell->rho > 0) {
            cell->v_avg /= cell->rho;
        }
      }
    }
  }
}



void calculate_unit_cost(Grid &grid) {
  /* iterate over each of 4 directions of each cell
     and calculate f_{M->i} and {C_M->i} where i is in {N, E, S, W} */

  for (auto &row : grid.grid) {
    for (auto &cell : row) {
      for (int dir = East; dir <= South; dir++) {
        /* dir = direction we're going TOWARDS */
        Edge *e = cell.edges[dir];

        /* Equation 8: topographical speed */
        float slope = ((e->h_grad * n_theta_int[dir]) - grid.s_min) /
            (grid.s_max - grid.s_min);
        float topo_speed = grid.f_max + slope * (grid.f_min - grid.f_max);

        /* Equation 9: flow speed */
        float flow_speed;
        if (cell.neighbors[dir]) {
          flow_speed = glm::dot(cell.neighbors[dir]->v_avg, n_theta_vec[dir]);
          if (flow_speed < 0.0f) flow_speed = 0.0f;
        } else {
          flow_speed = 0.0f;
        }

        if (!cell.neighbors[dir]) {
          cell.f[dir] = 0.0f;
        } else if (cell.neighbors[dir]->rho <= grid.rho_min) {
          cell.f[dir] = topo_speed;
        } else if (cell.neighbors[dir]->rho >= grid.rho_max) {
          cell.f[dir] = flow_speed;
        } else {
          /* Equation 10: linearly interpolate */
          float speed = topo_speed + (cell.neighbors[dir]->rho - grid.rho_min) /
              (grid.rho_max - grid.rho_min) * (flow_speed - topo_speed);
          cell.f[dir] = speed;
        }

        /* Equation 4: Cost field */
        if (cell.neighbors[dir]) {
          float cost = (grid.alpha * cell.f[dir] + grid.beta + grid.gamma * cell.neighbors[dir]->g) / cell.f[dir];
          cell.C[dir] = cell.f[dir] > 0 ? cost : std::numeric_limits<float>::infinity();
        } else {
          cell.C[dir] = std::numeric_limits<float>::infinity();
        }
      }
    }
  }
}

std::vector<Cell*> flatten(Grid &grid) {
  std::vector<Cell*> result;
  for (int j = 0; j < grid.height; j++) {
    for (int i = 0; i < grid.width; i++) {
      result.push_back(grid.getCell(i, j));
    }
  }
  return result;
}


void finite_differences_approx(Cell &cell) {
  float phi_mx;
  float phi_my;
  float mx;
  float my;
  Direction d_mx;
  Direction d_my;

  //Check boundary cases
  if (cell.neighbors[West] == nullptr) {
    phi_mx = cell.neighbors[East]->phi;
    mx = phi_mx + cell.C[East];
    d_mx = East;
  } else if (cell.neighbors[East] == nullptr) {
    phi_mx = cell.neighbors[West]->phi;
    mx = phi_mx + cell.C[West];
    d_mx = West;
  } else {
    //Otherwise, choose minimum phi + index between west/east directions
    float wx = cell.neighbors[West]->phi + cell.C[West];
    float ex = cell.neighbors[East]->phi + cell.C[East];

    if (wx < ex) {
      phi_mx = cell.neighbors[West]->phi;
      mx = phi_mx + cell.C[West];
      d_mx = West;
    } else {
      phi_mx = cell.neighbors[East]->phi;
      mx = phi_mx + cell.C[East];
      d_mx = East;
    }
  }
  //Check boundary cases
  if (cell.neighbors[South] == nullptr) {
    phi_my = cell.neighbors[North]->phi;
    my = phi_my + cell.C[North];
    d_my = North;
  } else if (cell.neighbors[North] == nullptr) {
    phi_my = cell.neighbors[South]->phi;
    my = phi_my + cell.C[South];
    d_my = South;
  } else {
    //Otherwise, choose minimum phi + index between north/south directions
    float ny = cell.neighbors[North]->phi + cell.C[North];
    float sy = cell.neighbors[South]->phi + cell.C[South];

    if (sy < ny) {
      phi_my = cell.neighbors[South]->phi;
      my = phi_my + cell.C[South];
      d_my = South;
    } else {
      phi_my = cell.neighbors[North]->phi;
      my = phi_my + cell.C[North];
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
    cell.edges[d_my]->phi_grad = phi_my - phi_m;
    cell.edges[d_my]->v = cell.edges[d_my]->phi_grad *
        (float) -cell.neighbors[d_my]->f[(d_my + 2) % 4];
  } else if (std::isinf(phi_my)) {
    float det = c_my;
    phi_m = phi_mx + std::sqrt(det);
    cell.edges[d_mx]->phi_grad = phi_mx - phi_m;
    cell.edges[d_mx]->v = cell.edges[d_mx]->phi_grad *
        (float) -cell.neighbors[d_mx]->f[(d_mx + 2) % 4];
  } else {
    float a = (c_mx + c_my);
    float b = -2.0f * (c_my * phi_mx + c_mx * phi_my) ;
    phi_m = (-1.0f * b) / (2.0f * a);

    cell.edges[d_mx]->phi_grad = phi_mx - phi_m;
    cell.edges[d_mx]->v = cell.edges[d_mx]->phi_grad *
        (float) -cell.neighbors[d_mx]->f[(d_mx + 2) % 4];
    cell.edges[d_my]->phi_grad = phi_my - phi_m;
    cell.edges[d_my]->v = cell.edges[d_my]->phi_grad *
        (float) -cell.neighbors[d_my]->f[(d_my + 2) % 4];
  }

  cell.phi_tmp = phi_m;
}

// Compare function for the heap
bool cmp(const Cell * a, const Cell * b) {
  return a->phi_tmp > b->phi_tmp;
}


void construct_dynamic_potential_field(Grid &grid, Group &group) {
  /* use fast-marching algorithm to calculate the potentials
     and potential gradients (phi and del phi)
     also calculate the velocity field of velocities v */
  glm::ivec2 goal = group.goal;
  Cell *goal_cell = grid.getCell(goal);
  std::vector<Cell*> candidates = {goal_cell};
  std::make_heap(candidates.begin(), candidates.end(), cmp);
  goal_cell->phi_tmp = 0;
  for (int i = 0; i < grid.getHeight() * grid.getWidth(); i++) {
    std::pop_heap(candidates.begin(), candidates.end(), cmp);
    Cell *curr = candidates.back();
    candidates.pop_back();

    curr->phi = curr->phi_tmp;
    curr->status = KNOWN;

    for (Cell *c : curr->neighbors) {
      if (c != nullptr && c->status != KNOWN) {
        finite_differences_approx(*c);
        if (c->status == UNKNOWN) {
          c->status = CANDIDATE;
          candidates.push_back(c);
          std::push_heap(candidates.begin(), candidates.end(), cmp);
        }
      }
    }
  }

}

void calc_phi_grad(Grid &grid) {
  /* Calculate phi gradients in Edges from phi values in Cells */
  for (auto &row : grid.grid) {
    // Leftmost edges
    row[0].edges[West]->phi_grad = 0.0f;
    for (auto &cell : row) {
      cell.edges[East]->phi_grad = cell.neighbors[East] ?
          cell.neighbors[East]->phi - cell.phi : 0.0f;
      cell.edges[North]->phi_grad = cell.neighbors[North] ?
          cell.neighbors[North]->phi - cell.phi : 0.0f;
    }
  }
  // Bottom row
  for (auto &cell : grid.grid[0]) {
    cell.edges[South]->phi_grad = 0.0f;
  }
}

void normalize_gradients(Grid &grid) {
  /* Take each cell and make sure the phi_grads at its Edges are normalized
   * relative to each other */
  for (auto &row: grid.grid) {
    for (auto &cell : row) {
      float x_diff = cell.edges[East]->phi_grad - cell.edges[West]->phi_grad;
      float y_diff = cell.edges[North]->phi_grad - cell.edges[South]->phi_grad;

      glm::vec2 normalized = glm::normalize(glm::vec2(x_diff, y_diff));

      float x_mult = x_diff != 0.0f ? normalized.x / x_diff : 1.0f;
      float y_mult = y_diff != 0.0f ? normalized.y / y_diff : 1.0f;

      cell.edges[East]->phi_grad *= x_mult;
      cell.edges[West]->phi_grad *= x_mult;
      cell.edges[North]->phi_grad *= y_mult;
      cell.edges[South]->phi_grad *= y_mult;
    }
  }
}

glm::vec2 interpolateTwo(float x, float x1, float x2, glm::vec2 v1, glm::vec2 v2) {
  //  return v1 * (cell_index + 1.5f - personLoc) + v2 * (personLoc - cell_index + 0.5f);
  return v1 * (x2 - x) + v2 * (x - x1); // assume x2-x1 = 1
}

glm::vec2 interpolateFour(float x, float y, float x1, float y1, float x2, float y2,
                          glm::vec2 v11, glm::vec2 v12, glm::vec2 v21, glm::vec2 v22) {
  // assume x2-x1 = 1 and y2-y1 = 1
  float x_weight1 = x2 - x;
  float x_weight2 = x - x1;
  glm::vec2 interp1 = v11 * x_weight1 + v21 * x_weight2;
  glm::vec2 interp2 = v12 * x_weight1 + v22 * x_weight2;
  glm::vec2 interp_final = interp1 * (y2 - y) + interp2 * (y - y1);
  return interp_final;
}

void crowd_advection(Grid &grid, Group &group) {
  /* update each person's position by interpolating into the vector field */

  /* Put velocity in each Cell from its edges */
  for (auto &row : grid.grid) {
    for (auto &cell : row) {
      cell.v_avg = glm::vec2((cell.edges[East]->v + cell.edges[West]->v) / 2.0f,
                             (cell.edges[North]->v + cell.edges[South]->v) / 2.0f);
    }
  }

  /* Set people positions */
  glm::vec2 velocity;
  glm::vec2 cellAvel, cellBvel, cellCvel, cellDvel;
  std::vector<Person> &people = group.people;
  for (auto &person : people) {
    // bilinearly interpolate velocity for this person relative to center of cells
    cellAvel = cellBvel = cellCvel = cellDvel = glm::vec2(0.0f);
    Cell *cellA = person.getCell(grid);
    Cell *cellB = cellA->neighbors[East];
    Cell *cellC = cellB ? cellA->neighbors[East]->neighbors[North] : nullptr;
    Cell *cellD = cellA->neighbors[North];

    if (cellA) {
      cellAvel = cellA->v_avg;
    }
    if (cellB) {
      cellBvel = cellB->v_avg;
    }
    if (cellC) {
      cellCvel = cellC->v_avg;
    }
    if (cellD) {
      cellDvel = cellD->v_avg;
    }

    if (!cellA && !cellD) {
      // interpolate y-axis of B,C
      velocity = interpolateTwo(person.getPos().y, cellB->j + 0.5f, cellC->j + 0.5f, cellBvel, cellCvel);
      velocity.x = std::max(velocity.x, 0.0f);
    } else if (!cellC && !cellB) {
      // interpolate y-axis of A,D
      velocity = interpolateTwo(person.getPos().y, cellA->j + 0.5f, cellD->j + 0.5f, cellAvel, cellDvel);
      velocity.x = std::min(velocity.x, 0.0f);
    } else if (!cellA && !cellB) {
      // interpolate x-axis of D,C
      velocity = interpolateTwo(person.getPos().y, cellD->i + 0.5f, cellC->i + 0.5f, cellDvel, cellCvel);
      velocity.y = std::max(velocity.y, 0.0f);
    } else if (!cellC && !cellD) {
      // interpolate x-axis of A,B
      velocity = interpolateTwo(person.getPos().y, cellA->i + 0.5f, cellB->i + 0.5f, cellAvel, cellBvel);
      velocity.y = std::min(velocity.y, 0.0f);
    } else {
      // interpolate between all 4 cells
      velocity = interpolateFour(person.getPos().x, person.getPos().y,
                                 cellA->i + 0.5f, cellA->j + 0.5f,
                                 cellC->i + 0.5f, cellC->j + 0.5f,
                                 cellAvel, cellBvel, cellDvel, cellCvel);
    }

    // update person's velocity and position
    person.setVelocity(velocity);
    person.setPos(person.getPos() + velocity);

//    // set velocity with current cell's average velocity only
//    person.setVelocity(person.getCell(grid)->v_avg);
//    person.setPos(person.getPos() + person.getCell(grid)->v_avg);
  }
}



void enforce_minimum_distance(Grid &grid, std::vector<Group> &groups) {
  /* iterate over all pairs in a threshold distance and push people
     apart symmetrically until min. distance is reached. may instead use a
     neighbor grid instead of the vector of Persons. */

  // generate spatial map
  grid.build_neighbor_map(groups);

  for (Group &group : groups) {
    // run pair-wise minimum distance enforcement for each person
    std::vector<Person> &people = group.people;
    for (auto &person : people) {
      grid.handle_collisions(person);
    }
  }
}

void test_structures() {
  std::cout << "Hello World" << std::endl;

  Grid grid(4, 3);

  Person dalton(2.6f, 1.3f, 0, 0, -5.0f);
  Cell *daltonCell = dalton.getCell(grid); // 2,1
  daltonCell->g = 1234;
  daltonCell->edges[North]->v = 420.0f;
  Cell *daltonAbove = grid.getCell(dalton.getGridIndex() + glm::ivec2(0, 1));
  std::cout << daltonAbove->edges[South]->v << std::endl;
  std::cout << daltonAbove->neighbors[South]->g << std::endl;

  float smth = 7 * n_theta_int[South];
  std::cout << "-7: " << smth << std::endl;
}

void test_potential_field() {
  Grid grid(4, 3);

  Person dalton(2.6f, 1.3f, 0.0f, 0.0f, -5.0f);
  std::vector<Person> people = {dalton};
  Group group(glm::vec2(3.2, 2.6), people);
  Cell *daltonCell = dalton.getCell(grid); // 2,1
  daltonCell->g = 1234;
  daltonCell->edges[North]->v = 420.0f;
  Cell *daltonAbove = grid.getCell(dalton.getGridIndex() + glm::ivec2(0, 1));
  std::cout << daltonAbove->edges[South]->v << std::endl;
  std::cout << daltonAbove->neighbors[South]->g << std::endl;

  construct_dynamic_potential_field(grid, group);
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
  std::cout << setprecision(2);
  // test_structures();
  // test_potential_field();

  /* Init point display */
  pointDisplay(argc, argv);



  /* Write points to a file */
  //write_points(points, filename);

  if (argc < 2) {
    std::cerr << "Please specify a configuration file" << std::endl;
    return -1;
  }
  json j;

  if (load_config(j, argv[1])) {
    std::cerr << "Error loading config file" << std::endl;
    return -1;
  }

  Grid grid(j);
  std::vector<Group> groups;
  if (load_groups(j, &groups)) {
    std::cerr << "Error loading groups from config file" << std::endl;
    return -1;
  }

  int iterations = j["iterations"];
  std::vector<std::vector<glm::vec2>> point_traj;
  for (int i = 0; i < iterations; i++) {
    std::cout << "Iteration " << i << std::endl;
    density_conversion(grid, groups);
    calculate_unit_cost(grid);
    for (Group &group : groups) {
      grid.clearGridVals();
      construct_dynamic_potential_field(grid, group);
      crowd_advection(grid, group);
      for (Person &p : group.people) {
        std::cout << p.getPos()[0] << " " << p.getPos()[1] << std::endl;
      }
      point_traj.push_back(points_from_groups(groups));

    }
    enforce_minimum_distance(grid, groups);

    if (i == 30) {
      set_points(point_traj);
      /* Display points function, replace with actual point vector */
      //display_points(grid.getWidth(), grid.getHeight());
    }
  }
  set_points(point_traj);
  /* Display points function, replace with actual point vector */
  display_points(grid.getWidth(), grid.getHeight());


  return 0;
}
