#include <iomanip>

#include "NBodySimulation.h"
#include <omp.h>
/**
 * You can compile this file with
 *   make step-1-g++   // Uses the GNU Compiler Collection.
 *   make step-1-icpx  // Uses the Intel compiler.
 * and run it with
 *   ./setp-1-g++
 *   ./step-1-icpx
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

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

  // Code that initialises and runs the simulation.
  NBodySimulation nbs;
  double startTime = omp_get_wtime();
  nbs.setUp(argc, argv);
  double endTime = omp_get_wtime();         // Capture end time
  double elapsedTime = endTime - startTime; // Calculate elapsed time

  // Capture end time
  // Calculate elapsed time

  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();
  double startTime1 = omp_get_wtime();

  while (!nbs.hasReachedEnd())
  {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  double endTime1 = omp_get_wtime();
  double elapsedTime1 = endTime1 - startTime1;
  nbs.printSummary();
  nbs.closeParaviewVideoFile();
  // std::cout << "update body execution time: " << elapsedTime1 << " seconds" << std::endl;
  // std::cout << "setup execution time: " << elapsedTime << " seconds" << std::endl;

  return 0;
}
