#include <cmath>

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

class NBodySimulation {

 private:
  double t;
  double tFinal;
  double tPlot;
  double tPlotDelta;

  int NumberOfBodies;

  /**
   * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
   * each pointer represents one molecule/particle/body.
   */
  double** x;

  /**
   * Equivalent to x storing the velocities.
   */
  double** v;

  /**
   * One mass entry per molecule/particle.
   */
  double*  mass;

  /**
   * Global time step size used.
   */
  double timeStepSize;

  /**
   * Maximum velocity of all particles.
   */
  double maxV;

  /**
   * Minimum distance between two elements.
   */
  double minDx;

  /**
   * Stream for video output file.
   */
  std::ofstream videoFile;

  /**
   * Output counters.
   */
  int snapshotCounter;
  int timeStepCounter;


 public:
  NBodySimulation ();
  ~NBodySimulation ();

  /**
   * Check that the number command line parameters is correct.
   */
  void checkInput(int argc, char** argv);

  /**
   * Set up scenario from the command line.
   *
   * If you need additional helper data structures, you can initialise them
   * here. Alternatively, you can introduce a totally new function to initialise
   * additional data fields and call this new function from main after setUp().
   * Either way is fine.
   *
   * The semantics of this operations are not to be changed in the assignment.
   */
  void setUp (int argc, char** argv);

  /**
   * Compute forces.
   *
   * The current force is gravity, i.e. (x_i-x_j) * m_i * m_j/r^3.
   **/
  double force_calculation (int i, int j, int direction);

  /**
   * Implement timestepping scheme and force updates.
   */
  void updateBody ();

  /**
   * Check if the last time step has been reached (simulation is completed).
   *
   * This operation is not to be changed in the assignment.
   */
  bool hasReachedEnd ();

  /**
   * Take simulations snapshopts and print summary to standard output.
   *
   * This operation is not to be changed in the assignment.
   */
  void takeSnapshot ();

  /**
   * Handle Paraview output.
   *
   * These operations are not to be changed in the assignment.
   *
   * The file format is documented at
   * http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
   */
  void openParaviewVideoFile ();
  void closeParaviewVideoFile ();
  void printParaviewSnapshot ();

  /**
   * Handle terminal output.
   *
   * These operations are not to be changed in the assignment.
   */
  void printSnapshotSummary ();
  void printSummary ();

};
