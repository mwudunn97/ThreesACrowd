//
// Created by Dalton on 2/24/2019.
//

#include <glm/glm.hpp>
#include "Grid.h"

Cell::Cell(Edge *edgeE, Edge *edgeN, Edge *edgeW, Edge *edgeS) {
  edges[East] = edgeE;
  edges[North] = edgeN;
  edges[West] = edgeW;
  edges[South] = edgeS;
}

void Grid::fill() {
  /* Create edges first */
  for (int i = 0; i < 2 * width * height + width + height; i++) {
    edges.emplace_back();
  }

  int edgeRowLen = 2 * width + 1;

  /* Create cells and connect to existing edges */
  for (int j = 0; j < width; j++) {
    grid.emplace_back();
    for (int i = 0; i < height; i++) {
      int edgeS_idx = j * edgeRowLen + i;
      grid[j].emplace_back(
          &edges[edgeS_idx + width + 1],
          &edges[edgeS_idx + edgeRowLen],
          &edges[edgeS_idx + width],
          &edges[edgeS_idx]);
    }
  }

  /* Attach neighbors */
  for (int j = 0; j < width; j++) {
    for (int i = 0; i < height; i++) {
      if (i > 0) grid[j][i].neighbors[West] = &grid[j][i-1];
      if (i < width - 1) grid[j][i].neighbors[East] = &grid[j][i+1];
      if (j < height - 1) grid[j][i].neighbors[North] = &grid[j+1][i];
      if (j > 0) grid[j][i].neighbors[South] = &grid[j-1][i];
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
  return x * 31 + y;
}

// construct a spatial map of neighbors for all people
void Grid::build_neighbor_map(std::vector<Person> &people) {
  // clear entries
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // hash each person into the map
  for (auto &person : people) {
    float hash = hash_position(person.getPos());
    if (map[hash] == NULL) {
      map[hash] = new std::vector<Person *>();
    }
    map[hash]->push_back(&person);
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
        total_correction += (1.0 - glm::length(dir)) * (dir / glm::length(dir));
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
  f_min = j["f_min"];
  f_max = j["f_max"];

  fill();
}

int Grid::getWidth() const {
  return width;
}

int Grid::getHeight() const {
  return height;
}

Cell *Grid::getCell(int i, int j) {
  return &grid[j][i];
}

Cell* Grid::getCell(glm::ivec2 ij) {
  return &grid[ij[1]][ij[0]];
}
