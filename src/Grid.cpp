//
// Created by Dalton on 2/24/2019.
//

#include "Grid.h"

Cell::Cell(Edge *edgeE, Edge *edgeN, Edge *edgeW, Edge *edgeS) {
    edges[East] = edgeE;
    edges[North] = edgeN;
    edges[West] = edgeW;
    edges[South] = edgeS;
}

Grid::Grid(int width, int height)
: width(width), height(height) {
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