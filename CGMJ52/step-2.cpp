#include <iomanip>
#include <iostream>
#include "NBodySimulation.h"

/**
 * You can compile this file with
 *   make step-2-g++   // Uses the GNU Compiler Collection.
 *   make step-2-icpx  // Uses the Intel compiler.
 * and run it with
 *   ./step-2-g++
 *   ./step-2-icpx
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

class BarnesHuttNBS : public NBodySimulation
{
public:
  struct vec
  {
    double x, y, z;

    // Constructor to initialize the vector
    vec(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}

    // Operator overloads for vector addition and scalar division (useful for force calculations)
    vec &operator+=(const vec &rhs)
    {
      x += rhs.x;
      y += rhs.y;
      z += rhs.z;
      return *this;
    }

    vec operator/(double scalar) const
    {
      return {x / scalar, y / scalar, z / scalar};
    }
  };

  vec *forces = new vec[NumberOfBodies];

  struct BHTreeNode
  {
    double mass = 0;
    vec position = {0, 0, 0};
    vec COM = {0, 0, 0};
    BHTreeNode *children[8];
    double size = 0;
    bool IsLeaf = true;
    int body = -1; // store the index of the body if leaf
  };

  BHTreeNode *root;

  BHTreeNode *BuildTree()
  {
    // create the bounds for x,y,z
    double minX = x[0][0], maxX = x[0][0], minY = x[0][1], maxY = x[0][1], minZ = x[0][2], maxZ = x[0][2], diameter;
    for (int i = 1; i < NumberOfBodies; i++)
    {
      if (x[i][0] < minX)
        minX = x[i][0];
      if (x[i][1] < minY)
        minY = x[i][1];
      if (x[i][2] < minZ)
        minZ = x[i][2];
      if (x[i][0] > maxX)
        maxX = x[i][0];
      if (x[i][1] > maxY)
        maxY = x[i][1];
      if (x[i][2] > maxZ)
        maxZ = x[i][2];
    }
    // create the largest size to make our quadrants from
    diameter = maxX - minX;
    if ((maxY - minY) > diameter)
      diameter = maxY - minY;
    if ((maxZ - minZ) > diameter)
      diameter = maxY - minY;

    root = new BHTreeNode;
    root->position.x = (maxX + minX) / 2.0;
    root->position.y = (maxY + minY) / 2.0;
    root->position.z = (maxZ + minZ) / 2.0;
    root->size = diameter;
    // root setup complete now need to add children
    for (int i = 0; i < NumberOfBodies; i++)
    {
      InsertToNode(i, root);
    }
    // reset the tree/ make a root node
    //  for all particles insert them into the tree
    return root;
  };
  // function to give the space which a body belongs to in the current tree
  char GetOct(vec TreePosition, int body)
  {
    char oct = 0;

    if (x[body][0] > TreePosition.x)
    {
      oct |= 1; // Bit 0 for x-axis
    }
    if (x[body][1] > TreePosition.y)
    {
      oct |= 2; // Bit 1 for y-axis
    }
    if (x[body][2] > TreePosition.z)
    {
      oct |= 4; // Bit 2 for z-axis
    }
    return oct;
  }

  void createChild(int child, BHTreeNode *node)
  {
    double childSize = node->size / 2.0;
    vec childCenter;
    childCenter.x = node->position.x + childSize * ((child & 1) ? 0.5 : -0.5);
    childCenter.y = node->position.y + childSize * ((child & 2) ? 0.5 : -0.5);
    childCenter.z = node->position.z + childSize * ((child & 4) ? 0.5 : -0.5);

    if (node->position.x == childCenter.x || node->position.y == childCenter.y || node->position.y == childCenter.y)
    {
      std::cout << "node and child have the same positions" << childCenter.x << ", "
                << childCenter.y << ", "
                << childCenter.z << ") with size: "
                << childSize << node->position.x << ", "
                << node->position.y << ", "
                << node->position.z << ") with size: "
                << node->size << std::endl;
    }
    node->children[child] = new BHTreeNode;
    node->children[child]->position = childCenter;
    node->children[child]->size = childSize;
    node->children[child]->IsLeaf = true;
  }

  void InsertToNode(int newBody, BHTreeNode *node)
  {
    if (node->IsLeaf)
    {
      if (node->body == -1)
      {

        // Leaf node without a body, directly store the new body here
        node->body = newBody;
      }
      else
      {

        // Leaf node already contains a body, subdivide and reinsert both bodies

        int oldChildIdx = GetOct(node->position, node->body);
        int newChildIdx = GetOct(node->position, newBody);
        createChild(newChildIdx, node);
        createChild(oldChildIdx, node);
        InsertToNode(node->body, node->children[oldChildIdx]);
        InsertToNode(newBody, node->children[newChildIdx]);

        node->IsLeaf = false; // The node is no longer a leaf after subdivision
        node->body = -1;      // Reset body index as it now contains more than one body
      }
    }
    else
    {

      // Non-leaf node, find the correct child for the new body and insert
      int childIdx = GetOct(node->position, newBody);
      if (node->children[childIdx] == nullptr)
      {
        createChild(childIdx, node); // Proper initialization needed
      }

      InsertToNode(newBody, node->children[childIdx]);
    }

    // After insertion, update this node's mass and COM

    // If the node contains a body or has children, calculate the new mass and COM
    node->mass += mass[newBody]; // Add the new body's mass
    // Recalculate COM

    double totalX = node->COM.x * node->mass + x[newBody][0] * mass[newBody];
    double totalY = node->COM.y * node->mass + x[newBody][1] * mass[newBody];
    double totalZ = node->COM.z * node->mass + x[newBody][2] * mass[newBody];

    node->COM.x = totalX / node->mass;
    node->COM.y = totalY / node->mass;
    node->COM.z = totalZ / node->mass;
  }

  vec CalculateForce(int body, double **x, BHTreeNode *node, float theta)
  {
    //////////////ended here unable to access x propperly

    std::cout << "x[0]: [" << x[0] << ", " << x[0] << ", " << x[0] << "]" << std::endl;

    vec force = {0, 0, 0};
    // node has no body
    if (node == nullptr || node->mass == 0 || (node->IsLeaf && node->body == -1))
    {
      std::cout << "node is empty or has no body in it or no mass " << body << std::endl;
      return force;
    } // node is a leaf and contains a body that is different to current body
    if (node->IsLeaf && node->body != -1 && node->body != body)
    {
      std::cout << "node is a leaf and is a different body " << body << std::endl;
      force = {force_calculation(body, node->body, 0),
               force_calculation(body, node->body, 1),
               force_calculation(body, node->body, 2)};
      return force;
    }
    std::cout << "starting distance calc " << body << std::endl;
    const double distance = 3;
    /*
    std::cout << "x[0]: [" << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << "]" << std::endl;

    const double distance = sqrt(
        (x[body][0] - node->COM.x) * (x[body][0] - node->COM.x) +
        (x[body][1] - node->COM.y) * (x[body][1] - node->COM.y) +
        (x[body][2] - node->COM.z) * (x[body][2] - node->COM.z));
        */
    if (node->size / distance < theta)
    {
      std::cout << "node is far enough away to be treated as point mass " << body << std::endl;
      double distance3 = distance * distance * distance + 1e-10; // Add small value to avoid division by zero
      force = {(x[body][0] - node->COM.x) * mass[body] * node->mass / distance3,
               (x[body][1] - node->COM.y) * mass[body] * node->mass / distance3,
               (x[body][2] - node->COM.z) * mass[body] * node->mass / distance3};
    }
    else // node is not far enough away to be treated as point mass so find force of individual parts
    {
      for (int i = 0; i < 8; i++)
      {
        if (node->children[i] != nullptr)
        {
          std::cout << "initiating recursion into children " << body << std::endl;
          force += CalculateForce(body, x, node->children[i], theta);
        }
        else
        {
          std::cout << "line is wrong here" << std::endl;
        }
      }
    }
    return force;
  }

  void BHTUpdateBody(float theta)
  {
    timeStepCounter++;
    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();
    BHTreeNode *root = BuildTree();
    for (int j = 0; j < NumberOfBodies; j++)
    {

      forces[j] = vec();
      forces[j] = CalculateForce(j, x, root, theta);
    }
    std::cout << "finnished force calc" << std::endl;
    // seperate blocks so that we are not changing positions of bodies
    // while still calculating interactions with other bodies
    for (int j = 0; j < NumberOfBodies; j++)
    {
      std::cout << "v[0] " << v[0] << std::endl;
      std::cout << "x[0] " << x[0] << std::endl;
      x[j][0] = x[j][0] + timeStepSize * v[j][0];
      x[j][1] = x[j][1] + timeStepSize * v[j][1];
      x[j][2] = x[j][2] + timeStepSize * v[j][2];

      v[j][0] = v[j][0] + timeStepSize * forces[j].x / mass[j];
      v[j][1] = v[j][1] + timeStepSize * forces[j].y / mass[j];
      v[j][2] = v[j][2] + timeStepSize * forces[j].z / mass[j];
      double V = std::sqrt(v[j][0] * v[j][0] + v[j][1] * v[j][1] + v[j][2] * v[j][2]);
      if (maxV < V)
      {
        maxV = V;
      }
    }
  };
};
/**
 * Main routine.
 *
 * No major changes are needed in the assignment. You can add initialisation or
 * or remove input checking, if you feel the need to do so. But keep in mind
 * that you may not alter what the program writes to the standard output.
 */

int main(int argc, char **argv)
{

  std::cout << std::setprecision(15);
  float theta = 1.0;
  // Code that initialises and runs the simulation.
  BarnesHuttNBS nbs;
  nbs.setUp(argc, argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();
  double testArray[1][3] = {{0.1, 0.2, 0.3}};
  nbs.BHTUpdateBody(theta);

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
  while (!nbs.hasReachedEnd())
  {
    nbs.BHTUpdateBody(theta);
    nbs.takeSnapshot();
  }
}
