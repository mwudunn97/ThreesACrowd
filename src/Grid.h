//
// Created by Dalton on 2/24/2019.
//

#ifndef THREESACROWD_GRID_H
#define THREESACROWD_GRID_H

#include <glm/vec4.hpp>

class Grid {
public:
    Grid();

    glm::vec4 getPoint();

private:
    glm::vec4 point;
};


#endif //THREESACROWD_GRID_H
