//
// Created by Dalton on 2/27/2019.
//

#ifndef THREESACROWD_PERSON_H
#define THREESACROWD_PERSON_H

#include <glm/vec2.hpp>
#include "Grid.h"

class Person {
public:
    Person();
    Person(float x, float y, float smelliness);

    /* Grid index of current position */
    glm::ivec2 getGridIndex() const;

    /* Cell of current position */
    Cell *getCell(Grid &grid);

private:
    float x;
    float y;
    float smelliness;
};


#endif //THREESACROWD_PERSON_H
