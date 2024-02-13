#include <iomanip>

#include "NBodySimulationVectorised.cpp"

/**
 * You can compile this file with
 *   make step-4-g++   // Uses the GNU Compiler Collection.
 *   make step-4-icpx  // Uses the Intel compiler.
 * and run it with
 *   ./step-4-g++
 *   ./step-4-icpx
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

class NBodySimulationParallelised : public NBodySimulationVectorised
{
#include "NBodySimulation.h"
#include <mm_malloc.h>
#include <omp.h>

public:
  int NumberOfBodies;
  int C;
  double *distance_x;
  double *distance_y;
  double *distance_z;
  double *velocity_x;
  double *velocity_y;
  double *velocity_z;
  double *force_x;
  double *force_y;
  double *force_z;
  double *mass;

  void checkInput(int argc, char **argv)
  {
    if (argc == 1)
    {
      std::cerr << "usage: " << std::string(argv[0])
                << " plot-time final-time dt objects" << std::endl
                << " Details:" << std::endl
                << " ----------------------------------" << std::endl
                << "  plot-time:       interval after how many time units to plot."
                   " Use 0 to switch off plotting"
                << std::endl
                << "  final-time:      simulated time (greater 0)" << std::endl
                << "  dt:              time step size (greater 0)" << std::endl
                << "  objects:         any number of bodies, specified by position, velocity, mass" << std::endl
                << std::endl
                << "Examples of arguments:" << std::endl
                << "+ One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
                << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0" << std::endl
                << "+ One body spiralling around the other" << std::endl
                << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0" << std::endl
                << "+ Three-body setup from first lecture" << std::endl
                << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0" << std::endl
                << "+ Five-body setup" << std::endl
                << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0" << std::endl
                << std::endl;

      throw -1;
    }
    else if ((argc - 4) % 7 != 0)
    {
      std::cerr << "error in arguments: each body is given by seven entries"
                   " (position, velocity, mass)"
                << std::endl;
      std::cerr << "got " << argc << " arguments"
                                     " (three of them are reserved)"
                << std::endl;
      std::cerr << "run without arguments for usage instruction" << std::endl;
      throw -2;
    }
  }

  void setUp(int argc, char **argv)
  {
    checkInput(argc, argv);

    NumberOfBodies = (argc - 4) / 7;
    /*
            mass = new double[NumberOfBodies];
            velocity_x = new double[NumberOfBodies];
            velocity_y = new double[NumberOfBodies];
            velocity_z = new double[NumberOfBodies];
            distance_x = new double[NumberOfBodies];
            distance_y = new double[NumberOfBodies];
            distance_z = new double[NumberOfBodies];
            force_x = new double[NumberOfBodies];
            force_y = new double[NumberOfBodies];
            force_z = new double[NumberOfBodies];
            */
    const int alignment = 64; // Common cache line size

    // Allocate memory for each property
    mass = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);
    velocity_x = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);
    velocity_y = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);
    velocity_z = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);
    distance_x = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);
    distance_y = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);
    distance_z = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);
    force_x = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);
    force_y = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);
    force_z = (double *)_mm_malloc(NumberOfBodies * sizeof(double), alignment);

    // Check for successful allocation
    if (!mass || !velocity_x || !velocity_y || !velocity_z || !distance_x || !distance_y || !distance_z || !force_x || !force_y || !force_z)
    {
      std::cout << "MM_Malloc not worked" << std::endl;
      exit(-1);
    }

    C = 0.01 / NumberOfBodies;
    int readArgument = 1;

    tPlotDelta = std::stof(argv[readArgument]);
    readArgument++;
    tFinal = std::stof(argv[readArgument]);
    readArgument++;
    timeStepSize = std::stof(argv[readArgument]);
    readArgument++;

    for (int i = 0; i < NumberOfBodies; i++)
    {

      distance_x[i] = std::stof(argv[readArgument]);
      readArgument++;
      distance_y[i] = std::stof(argv[readArgument]);
      readArgument++;
      distance_z[i] = std::stof(argv[readArgument]);
      readArgument++;

      velocity_x[i] = std::stof(argv[readArgument]);
      readArgument++;
      velocity_y[i] = std::stof(argv[readArgument]);
      readArgument++;
      velocity_z[i] = std::stof(argv[readArgument]);
      readArgument++;

      mass[i] = std::stof(argv[readArgument]);
      readArgument++;
    }
    for (int i = 0; i < NumberOfBodies; i++)
    {
      if (mass[i] <= 0.0)
      {
        std::cerr << "invalid mass for body " << i << std::endl;
        exit(-2);
      }
    }
    std::cout << "created setup with " << NumberOfBodies << " bodies"
              << std::endl;

    if (tPlotDelta <= 0.0)
    {
      std::cout << "plotting switched off" << std::endl;
      tPlot = tFinal + 1.0;
    }
    else
    {
      std::cout << "plot initial setup plus every " << tPlotDelta
                << " time units" << std::endl;
      tPlot = 0.0;
    }
  }

  void check_collision()
  {
    int *merge = new int[NumberOfBodies];
    std::fill_n(merge, NumberOfBodies, -1);

    int newNBodies = NumberOfBodies;
#pragma omp target parallel for
    for (int i = 0; i < NumberOfBodies - 1; i++)
    {

      for (int j = i + 1; j < NumberOfBodies; j++)
      {
        double distance = sqrt(
            (distance_x[j] - distance_x[i]) * (distance_x[j] - distance_x[i]) +
            (distance_y[j] - distance_y[i]) * (distance_y[j] - distance_y[i]) +
            (distance_z[j] - distance_z[i]) * (distance_z[j] - distance_z[i]));

        if (distance < C * (mass[i] + mass[j]))
        {
          // Simplified merge logic for clarity
          merge[j] = (merge[i] == -1) ? i : merge[i];
          newNBodies -= (merge[i] == -1) ? 1 : 0;
        }
      }
    }

    // Allocate new SoA arrays for position, velocity, and mass
    double *newDistance_x = new double[newNBodies];
    double *newDistance_y = new double[newNBodies];
    double *newDistance_z = new double[newNBodies];
    double *newVelocity_x = new double[newNBodies];
    double *newVelocity_y = new double[newNBodies];
    double *newVelocity_z = new double[newNBodies];
    double *newM = new double[newNBodies];

    int newIndex = 0;

    for (int i = 0; i < NumberOfBodies; i++)
    {
      if (merge[i] != -1)
        continue; // Skip bodies that are sources in a merge

      // Initialize mass and position/velocity with the body's current state
      newM[newIndex] = mass[i];
      newDistance_x[newIndex] = distance_x[i];
      newDistance_y[newIndex] = distance_y[i];
      newDistance_z[newIndex] = distance_z[i];
      newVelocity_x[newIndex] = velocity_x[i];
      newVelocity_y[newIndex] = velocity_y[i];
      newVelocity_z[newIndex] = velocity_z[i];

      // Apply merges

      for (int j = 0; j < NumberOfBodies; j++)
      {
        if (merge[j] == i)
        { // Body j merges into body i
          double totalMass = newM[newIndex] + mass[j];
          newDistance_x[newIndex] = (newDistance_x[newIndex] * newM[newIndex] + distance_x[j] * mass[j]) / totalMass;
          newDistance_y[newIndex] = (newDistance_y[newIndex] * newM[newIndex] + distance_y[j] * mass[j]) / totalMass;
          newDistance_z[newIndex] = (newDistance_z[newIndex] * newM[newIndex] + distance_z[j] * mass[j]) / totalMass;
          newVelocity_x[newIndex] = (newVelocity_x[newIndex] * newM[newIndex] + velocity_x[j] * mass[j]) / totalMass;
          newVelocity_y[newIndex] = (newVelocity_y[newIndex] * newM[newIndex] + velocity_y[j] * mass[j]) / totalMass;
          newVelocity_z[newIndex] = (newVelocity_z[newIndex] * newM[newIndex] + velocity_z[j] * mass[j]) / totalMass;
          newM[newIndex] = totalMass;
        }
      }

      newIndex++;
    }

    // Free old arrays
    delete[] distance_x;
    delete[] distance_y;
    delete[] distance_z;
    delete[] velocity_x;
    delete[] velocity_y;
    delete[] velocity_z;
    delete[] mass;

    // Update pointers to new arrays
    distance_x = newDistance_x;
    distance_y = newDistance_y;
    distance_z = newDistance_z;
    velocity_x = newVelocity_x;
    velocity_y = newVelocity_y;
    velocity_z = newVelocity_z;
    mass = newM;
    NumberOfBodies = newNBodies;

    delete[] merge;
  }
  void updateBody()

  {
    int i = 0;
    const int tileSize = 8;
    const int alignment = 64;
    timeStepCounter++;
    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();
#pragma omp target parallel for shared(distance_x, distance_y, distance_z, mass)
    for (int chunk = 0; chunk < NumberOfBodies; chunk += tileSize)
    {
      double f_xchunk[tileSize];
      double f_ychunk[tileSize];
      double f_zchunk[tileSize];
#pragma omp simd
      for (int j = 0; j < tileSize; j++)
      {
        f_xchunk[j] = 0.0;
        f_ychunk[j] = 0.0;
        f_zchunk[j] = 0.0;
      }

      for (int j = 0; j < NumberOfBodies; j++)
      {
        for (int i = chunk; i < chunk + tileSize; i++)
        {
          if (i == j)
            continue;
          double dx, dy, dz;
          double distance = 0.0;
          double distanceInv = 0.0;
          dx = distance_x[j] - distance_x[i];
          dy = distance_y[j] - distance_y[i];
          dz = distance_z[j] - distance_z[i];
          distance = sqrtl(dx * dx + dy * dy + dz * dz);
          distanceInv = 1 / distance;
          minDx = std::min(minDx, distance);

          f_xchunk[i - chunk] += dx * mass[j] * distanceInv * distanceInv * distanceInv;
          f_ychunk[i - chunk] += dy * mass[j] * distanceInv * distanceInv * distanceInv;
          f_zchunk[i - chunk] += dz * mass[j] * distanceInv * distanceInv * distanceInv;
        }
      }
#pragma omp simd
      for (int k = 0; k < tileSize; k++)
      {
        force_x[k + chunk] = f_xchunk[k];
        force_y[k + chunk] = f_ychunk[k];
        force_z[k + chunk] = f_zchunk[k];
      }
    }
#pragma omp target parallel for
    for (int j = 0; j < NumberOfBodies; j++)
    {

      distance_x[j] = distance_x[j] + timeStepSize * velocity_x[j];
      distance_y[j] = distance_y[j] + timeStepSize * velocity_y[j];
      distance_z[j] = distance_z[j] + timeStepSize * velocity_z[j];
      double massInv = mass[j];
      velocity_x[j] = velocity_x[j] + timeStepSize * force_x[j];
      velocity_y[j] = velocity_y[j] + timeStepSize * force_y[j];
      velocity_z[j] = velocity_z[j] + timeStepSize * force_z[j];
      double V = std::sqrt(velocity_x[j] * velocity_x[j] + velocity_y[j] * velocity_y[j] + velocity_z[j] * velocity_z[j]);
      if (maxV < V)
      {
        maxV = V;
      }
    }

    check_collision();

    t += timeStepSize;
  }

  bool
  hasReachedEnd()
  {
    return t > tFinal;
  }

  void takeSnapshot()
  {
    if (t >= tPlot)
    {
      printParaviewSnapshot();
      printSnapshotSummary();
      tPlot += tPlotDelta;
    }
  }

  void printParaviewSnapshot()
  {
    static int counter = -1;
    counter++;
    std::stringstream filename, filename_nofolder;
    filename << "paraview-output/result-" << counter << ".vtp";
    filename_nofolder << "result-" << counter << ".vtp";
    std::ofstream out(filename.str().c_str());
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float64\""
           " NumberOfComponents=\"3\""
           " format=\"ascii\">";

    for (int i = 0; i < NumberOfBodies; i++)
    {
      out << distance_x[i]
          << " "
          << distance_y[i]
          << " "
          << distance_z[i]
          << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>" << std::endl;

    out.close();

    videoFile << "<DataSet timestep=\"" << counter
              << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
              << "\"/>" << std::endl;
  }

  void printSummary()
  {
    std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
    std::cout << "Position of first remaining object: "
              << distance_x[0] << ", " << distance_y[0] << ", " << distance_z[0] << std::endl;
  }
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

  // Code that initialises and runs the simulation.
  NBodySimulationParallelised nbs;
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
  std::cout << "update body execution time: " << elapsedTime1 << " seconds" << std::endl;
  std::cout << "setup execution time: " << elapsedTime << " seconds" << std::endl;

  return 0;
}
