//
// Created by Dalton on 2/24/2019.
//

#ifndef THREESACROWD_GRID_H
#define THREESACROWD_GRID_H

#include <glm/vec4.hpp>
#include <glm/vec2.hpp>
#include <vector>
#include <array>

struct Edge {
    glm::vec2 g_grad;
    glm::vec2 phi_grad;
    glm::vec2 v;
};

struct Cell {
    float g;
    float phi;
    float rho;
    float h;
    glm::vec2 v_avg;
    /* Going INTO this cell: E,N,W,S */
    glm::vec4 f;
    glm::vec4 C;
    std::array<Edge*, 4> edges;
};

class Grid {
public:
    Grid();

private:
    std::vector<std::vector<Cell>> grid;
    std::vector<Edge> edges;
};


#endif //THREESACROWD_GRID_H
