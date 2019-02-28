//
// Created by Dalton on 2/27/2019.
//

#include "Person.h"

Person::Person()
: x(0.0f), y(0.0f), smelliness(0.0f) {}

Person::Person(float x, float y, float smelliness)
: x(x), y(y), smelliness(smelliness) {}

glm::ivec2 Person::getGridIndex() const {
    return glm::ivec2(static_cast<int>(x), static_cast<int>(y));
}

Cell *Person::getCell(Grid &grid) {
    return grid.getCell(static_cast<int>(x), static_cast<int>(y));
}