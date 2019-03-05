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

using json = nlohmann::json;

void density_conversion(Grid &grid, std::vector<Person> &people, double lambda) {
  // TODO 4.1: convert positions of Persons into densities and
  // insert into Grid. Also calculate average velocities of each cell.

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
  return;
}


void construct_dynamic_potential_field(Grid &grid) {
  // TODO 4.3: use fast-marching algorithm to calculate the potentials
  // and potential gradients (phi and del phi)
  // also calculate the velocity field of velocities v
  return;
}


void crowd_avection(Grid &grid, std::vector<Person> &people) {
  // TODO 4.4: update each person's position by interpolating into the vector
  // field
  for (auto &person : people) {
    person.setPos(person.getPos() + person.getCell(grid)->v_avg);
  }
}


void enforce_minimum_distsance(Grid &grid, std::vector<Person> &people) {
  // TODO 4.5: iterate over all pairs in a threshold distance and push people
  // apart symmetrically until min. distance is reached. may instead use a
  // neighbor grid instead of the vector of Persons.
  return;
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

int load_config(json &j, char *config) {
    std::ifstream f(config);
    if (!f.good()) {
        std::cerr << "Config file " << config << " does not exist" << std::endl;
        return -1;
    }
    f >> j;
    return 0;
}

int main(int argc, char* argv[]) {
    // test_structures();

    if (argc < 2) {
        std::cerr << "Please specify a configuration file" << std::endl;
        return -1;
    }
    json j;
    if (load_config(j, argv[1])) {
        return -1;
    }

    Grid grid(j["grid"]["width"], j["grid"]["height"]);
    std::cout << grid.getWidth() << std::endl;
}