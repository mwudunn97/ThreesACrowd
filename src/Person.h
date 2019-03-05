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
    Person(glm::vec2 pos, glm::vec2 velocity, float smelliness);
    Person(float x, float y, float vx, float vy, float smelliness);

    /* Grid index of current position */
    glm::ivec2 getGridIndex() const;

    /* Cell of current position */
    Cell *getCell(Grid &grid);

    glm::vec2 getPos();
    glm::vec2 getVelocity();

private:
    glm::vec2 pos;
    glm::vec2 velocity;

    float smelliness;
};


#endif //THREESACROWD_PERSON_H
