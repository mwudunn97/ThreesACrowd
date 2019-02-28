//
// Created by Dalton on 2/24/2019.
//

#ifndef THREESACROWD_GRID_H
#define THREESACROWD_GRID_H

#include <glm/vec4.hpp>
#include <glm/vec2.hpp>
#include <vector>
#include <array>

enum Direction {
    East = 0,
    North = 1,
    West = 2,
    South = 3
};

struct Edge {
    glm::vec2 g_grad;
    glm::vec2 phi_grad;
    glm::vec2 v;
};

struct Cell {
    Cell(Edge *edgeE, Edge *edgeN, Edge *edgeW, Edge *edgeS);

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

/* Row-major grid
 * 0,0 in bottom left */
class Grid {
public:
    Grid(int width, int height);

    int getWidth() const;
    int getHeight() const;

    Cell *getCell(int i, int j);
    Cell *getCell(glm::ivec2 ij);


    std::vector<std::vector<Cell>> grid;
    std::vector<Edge> edges;

    int width;
    int height;
};


#endif //THREESACROWD_GRID_H
