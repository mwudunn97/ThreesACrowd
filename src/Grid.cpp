//
// Created by Dalton on 2/24/2019.
//

#include "Grid.h"

Grid::Grid() {
    this->point = glm::vec4(0.0f, 1.0f, 2.0f, 3.0f);
}

glm::vec4 Grid::getPoint() {
    return this->point;
}