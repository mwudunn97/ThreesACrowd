//
// Created by Dalton on 2/27/2019.
//

#include "Person.h"
#include "Grid.h"

Person::Person()
: pos(0, 0), velocity(0, 0), goal(0, 0), smelliness(0.0f) {}

Person::Person(glm::vec2 pos, glm::vec2 velocity, glm::vec2 goal, float smelliness)
: pos(pos), velocity(velocity), goal(goal), smelliness(smelliness) {}

Person::Person(float x, float y, float vx, float vy, float gx, float gy, float smelliness)
: pos(x, y), velocity(vx, vy), goal(gx, gy), smelliness(smelliness) {}

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

glm::vec2 Person::getGoal() {
	return goal;
}

void Person::setPos(glm::vec2 pos) {
  this->pos = pos;
}

void Person::setVelocity(glm::vec2 pos) {
  this->velocity = velocity;
}

Group::Group(glm::vec2 goal, std::vector<Person> people)
 : goal(goal), people(std::move(people)) {}
