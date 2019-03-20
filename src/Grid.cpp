//
// Created by Dalton on 2/24/2019.
//

#include <glm/glm.hpp>
#include "Grid.h"
#include <exception>
#include <iostream>
#include <iomanip>

Cell::Cell(Edge *edgeE, Edge *edgeN, Edge *edgeW, Edge *edgeS, int i, int j)
 : i(i), j(j) {
  edges[East] = edgeE;
  edges[North] = edgeN;
  edges[West] = edgeW;
  edges[South] = edgeS;
  this->phi = std::numeric_limits<float>::infinity();
}


void Grid::fill() {
  /* Create edges first */
  for (int i = 0; i < 2 * width * height + width + height; i++) {
    edges.emplace_back();
  }

  int edgeRowLen = 2 * width + 1;

  /* Create cells and connect to existing edges */
  for (int j = 0; j < height; j++) {
    grid.emplace_back();
    for (int i = 0; i < width; i++) {
      int edgeS_idx = j * edgeRowLen + i;
      grid[j].emplace_back(
          &edges[edgeS_idx + width + 1],
          &edges[edgeS_idx + edgeRowLen],
          &edges[edgeS_idx + width],
          &edges[edgeS_idx],
          i, j);
    }

  }

  /* Attach neighbors */
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      grid[j][i].neighbors[West] = i > 0 ? &grid[j][i-1] : nullptr;
      grid[j][i].neighbors[East] = i < width - 1 ?  &grid[j][i+1] : nullptr;
      grid[j][i].neighbors[North] = j < height - 1 ? &grid[j+1][i] : nullptr;
      grid[j][i].neighbors[South] = j > 0 ?  &grid[j-1][i] : nullptr;
    }
  }
}

void Grid::clearGridVals() {

  for (int j = 0; j < width; j++) {
    for (int i = 0; i < height; i++) {
      Cell * cell = getCell(i, j);
      cell->phi = std::numeric_limits<float>::infinity();
      cell->phi_tmp = std::numeric_limits<float>::infinity();
      cell->status = UNKNOWN;
      for (int dir = East; dir <= South; dir++) {
        Edge *e = cell->edges[dir];
        e->phi_grad = 0.0f;
        e->v = 0.0f;
      }
    }

  }
}


// compute a hash based on a position
float Grid::hash_position(glm::vec2 pos) {
  // TODO: find the right constant for this
  double subdiv_factor = 5.0;

  // calculate bounds of this rectangle
  double w = std::min(1.0, width / subdiv_factor);
  double h = std::min(1.0, height / subdiv_factor);

  // calculate floor
  double x = floor(pos[0] / w);
  double y = floor(pos[0] / h);
  return (float) (x * 31 + y);
}

// construct a spatial map of neighbors for all people
void Grid::build_neighbor_map(std::vector<Group> &group) {
  // clear entries
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // hash each person into the map
  for (auto &group : group) {
      std::vector<Person> &people = group.people;
      for (auto &person : people) {
          float hash = hash_position(person.getPos());
          if (map[hash] == NULL) {
              map[hash] = new std::vector<Person *>();
          }
          map[hash]->push_back(&person);
      }
  }
}

// check for collisions between people and correct for it
void Grid::handle_collisions(Person &person) {
  float hash = hash_position(person.getPos());
  std::vector<Person*>* neighbors = map[hash];

  glm::vec2 total_correction(0.0f, 0.0f);
  float num_corrections = 0;
  if (neighbors != NULL) {
    for (Person* other : *neighbors) {
      if (other == &person) continue;

      // if these two persons collide, move this one away from the other
      glm::vec2 dir = person.getPos() - other->getPos();
      if (glm::length(dir) < 1.0) {
        glm::length(dir);
        total_correction += (1.0f - glm::length(dir)) * (dir / glm::length(dir));
        num_corrections++;
      }
    }
  }
  if (num_corrections > 0) {
    person.setPos(person.getPos() + total_correction / num_corrections);
  }
}



Grid::Grid(int width, int height)
: width(width), height(height) {
  fill();
}

Grid::Grid(json &j) {
  width =   j["grid"]["width"];
  height =  j["grid"]["height"];
  alpha =   j["alpha"];
  beta =    j["beta"];
  gamma =   j["gamma"];
  rho_min = j["rho_min"];
  rho_max = j["rho_max"];
  f_min =   j["f_min"];
  f_max =   j["f_max"];
  s_min =   j["s_min"];
  s_max =   j["s_max"];
  lambda =  j["lambda"];

  fill();
}

int Grid::getWidth() const {
  return width;
}

int Grid::getHeight() const {
  return height;
}

Cell *Grid::getCell(int i, int j) {
  if (i < 0 || i >= width || j < 0 || j >= height) {
    std::cerr << i << " " << j << std::endl;
    throw std::runtime_error("Cell out of bounds");
  }
  return &grid[j][i];
}

Cell* Grid::getCell(glm::ivec2 ij) {
  if (ij[0] < 0 || ij[0] >= width || ij[1] < 0 || ij[1] >= height) {
    std::cerr << ij[0] << " " << ij[1] << std::endl;
    throw std::runtime_error("Cell out of bounds");
  }
  return &grid[ij[1]][ij[0]];
}

Cell* Grid::getCellFromPos(glm::vec2 pos) {
  return getCell(static_cast<int>(pos[0]), static_cast<int>(pos[1]));
}

void Grid::print_v_avg() {
  for (auto row = grid.rbegin(); row != grid.rend(); ++row) {
    for (auto &cell : *row) {
      std::cout << std::setw(6) << cell.v_avg[0] << " " << std::setw(6) << cell.v_avg[1] << " ||";
    }
    std::cout << std::endl;
  }
}

void Grid::print_density() {
  for (auto row = grid.rbegin(); row != grid.rend(); ++row) {
    for (auto &cell : *row) {
      std::cout << std::setw(5) << cell.rho << " ";
    }
    std::cout << std::endl;
  }
}

void Grid::print_f() {
  for (auto row = grid.rbegin(); row != grid.rend(); ++row) {
    for (auto &cell : *row) {
      std::cout << cell.f[0] << "|" << cell.f[1] << "|" << cell.f[2] <<
                "|" << cell.f[3] << " ";
    }
    std::cout << std::endl;
  }
}

void Grid::print_C() {
  for (auto row = grid.rbegin(); row != grid.rend(); ++row) {
    for (auto &cell : *row) {
      std::cout << cell.C[0] << "|" << cell.C[1] << "|" << cell.C[2] <<
                "|" << cell.C[3] << " ";
    }
    std::cout << std::endl;
  }
}


void Grid::print_phi() {
  for (auto row = grid.rbegin(); row != grid.rend(); ++row) {
    for (auto &cell : *row) {
      std::cout << std::setw(4) << cell.phi << " ";
    }

    std::cout << std::endl;
  }
}

void Grid::print_phi_grad() {
  std::cout << "   ";
  for (auto &cell : grid[height - 1]) {
    std::cout << std::setw(7) << cell.edges[North]->phi_grad << " ";
  }
  std::cout << std::endl;
  for (auto row = grid.rbegin(); row != grid.rend(); ++row) {
    std::cout << std::setw(7) <<  (*row)[0].edges[West]->phi_grad << " ";
    for (auto &cell : *row) {
      std::cout << std::setw(7) << cell.edges[East]->phi_grad << " ";
    }
    std::cout << std::endl << "   ";
    for (auto &cell : *row) {
      std::cout << std::setw(7) << cell.edges[South]->phi_grad << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void Grid::print_v() {
  std::cout << "   ";
  for (auto &cell : grid[height - 1]) {
    std::cout << std::setw(7) << cell.edges[North]->v<< " ";
  }
  std::cout << std::endl;
  for (auto row = grid.rbegin(); row != grid.rend(); ++row) {
    std::cout << std::setw(7) <<(*row)[0].edges[West]->v << " ";
    for (auto &cell : *row) {
      std::cout << std::setw(7) << cell.edges[East]->v << " ";
    }
    std::cout << std::endl << "   ";
    for (auto &cell : *row) {
      std::cout << std::setw(7) << cell.edges[South]->v << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}