/*
 * ThreesACrowd - implementation of "Continuum Crowds", Treuille, et al.
 * https://grail.cs.washington.edu/projects/crowd-flows/78-treuille.pdf
 * Dalton Omens, Marc WuDunn, Jessie Yang
*/

#include <iostream>
#include "Grid.h"

int main(int argc, char* argv[]) {
    std::cout << "Hello World" << std::endl;

    Grid grid;
}

void density_conversion(Grid &grid, std::vector<Person> &people) {
  // TODO 4.1: convert positions of Persons into densities and
  // insert into Grid. Also calculate average velocities of each cell.
  return;
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
  return;
}


void enforce_minimum_distsance(Grid &grid, std::vector<Person> &people) {
  // TODO 4.5: iterate over all pairs in a threshold distance and push people
  // apart symmetrically until min. distance is reached. may instead use a
  // neighbor grid instead of the vector of Persons.
  return;
}
