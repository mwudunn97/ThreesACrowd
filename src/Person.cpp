//
// Created by Dalton on 2/27/2019.
//

#include "Person.h"

Person::Person()
: pos(0, 0), velocity(0, 0), smelliness(0.0f) {}

Person::Person(glm::vec2 pos, glm::vec2 velocity, float smelliness)
: pos(pos), velocity(velocity), smelliness(smelliness) {}

Person::Person(float x, float y, float vx, float vy, float smelliness)
: pos(x, y), velocity(vx, vy), smelliness(smelliness) {}

glm::ivec2 Person::getGridIndex() const {
  return glm::ivec2(static_cast<int>(pos[0]), static_cast<int>(pos[1]));
}

Cell *Person::getCell(Grid &grid) {
  return grid.getCell(static_cast<int>(pos[0]), static_cast<int>(pos[1]));
}

glm::vec2 Person::getPos() {
  return pos;
}

glm::vec2 Person::getVelocity() {
  return velocity;
}

void Person::setPos(glm::vec2 pos) {
  this->pos = pos;
}

void Person::setVelocity(glm::vec2 pos) {
  this->velocity = velocity;
}

Group::Group(glm::vec2 goal, std::vector<Person> people)
 : goal(goal), people(std::move(people)) {}
