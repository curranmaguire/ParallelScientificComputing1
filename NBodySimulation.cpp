#include "NBodySimulation.h"

NBodySimulation::NBodySimulation() : t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
                                     x(nullptr), v(nullptr), mass(nullptr),
                                     timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
                                     snapshotCounter(0), timeStepCounter(0), force1(nullptr), force2(nullptr), force3(nullptr), C(0){};

NBodySimulation::~NBodySimulation()
{
  if (x != nullptr)
  {
    for (int i = 0; i < NumberOfBodies; i++)
      delete[] x[i];
    delete[] x;
  }
  if (v != nullptr)
  {
    for (int i = 0; i < NumberOfBodies; i++)
      delete[] v[i];
    delete[] v;
  }

  if (mass != nullptr)
  {
    delete[] mass;
  }
  if (force1 != nullptr)
  {
    delete[] force1;
  }
  if (force2 != nullptr)
  {
    delete[] force2;
  }
  if (force3 != nullptr)
  {
    delete[] force3;
  }
}

void NBodySimulation::checkInput(int argc, char **argv)
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

void NBodySimulation::setUp(int argc, char **argv)
{
  checkInput(argc, argv);

  NumberOfBodies = (argc - 4) / 7;

  x = new double *[NumberOfBodies];
  v = new double *[NumberOfBodies];
  mass = new double[NumberOfBodies];
  force1 = new double[NumberOfBodies];
  force2 = new double[NumberOfBodies];
  force3 = new double[NumberOfBodies];
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
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]);
    readArgument++;
    x[i][1] = std::stof(argv[readArgument]);
    readArgument++;
    x[i][2] = std::stof(argv[readArgument]);
    readArgument++;

    v[i][0] = std::stof(argv[readArgument]);
    readArgument++;
    v[i][1] = std::stof(argv[readArgument]);
    readArgument++;
    v[i][2] = std::stof(argv[readArgument]);
    readArgument++;

    mass[i] = std::stof(argv[readArgument]);
    readArgument++;

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

double NBodySimulation::force_calculation(int i, int j, int direction)
{
  // Euclidean distance
  const double distance = sqrt(
      (x[j][0] - x[i][0]) * (x[j][0] - x[i][0]) +
      (x[j][1] - x[i][1]) * (x[j][1] - x[i][1]) +
      (x[j][2] - x[i][2]) * (x[j][2] - x[i][2]));
  const double distance3 = distance * distance * distance;
  minDx = std::min(minDx, distance);

  return (x[i][direction] - x[j][direction]) * mass[i] * mass[j] / distance3;
}

void NBodySimulation::check_collision()
{
  // Assuming `C` and other relevant variables are defined correctly
  int *merge = new int[NumberOfBodies];
  std::fill_n(merge, NumberOfBodies, -1);

  int newNBodies = NumberOfBodies;
  for (int i = 0; i < NumberOfBodies - 1; i++)
  {
    for (int j = i + 1; j < NumberOfBodies; j++)
    {
      if (i == j)
        continue;
      double distance = sqrt((x[j][0] - x[i][0]) * (x[j][0] - x[i][0]) +
                             (x[j][1] - x[i][1]) * (x[j][1] - x[i][1]) +
                             (x[j][2] - x[i][2]) * (x[j][2] - x[i][2])); // Calculate as before;
      if (distance < C * (mass[i] + mass[j]))
      {
        if (merge[i] == -1 && merge[j] == -1)
        {               // Neither body has been merged yet
          merge[j] = i; // Merge j into i
          newNBodies--; // Decrease count of bodies only when a new merge is found
        }
        else
        {
          if (merge[i] != -1)
          {
            merge[j] = merge[i]; // merge the bodies into the original
          }
          else
          {
            if (merge[j] != -1)
            {
              merge[i] = merge[j];
            }
            else
            {
              if (merge[i] != -1 && merge[j] != -1)
              {
                // if both are already merging then
              }
            }
          }
        }
      }
    }
  }
  double **newX = new double *[newNBodies];
  double **newV = new double *[newNBodies];
  double *newM = new double[newNBodies];

  int newIndex = 0;
  for (int i = 0; i < NumberOfBodies; i++)
  {
    // Skip bodies that are sources in a merge
    if (merge[i] != -1)
      continue;

    // Allocate new arrays for position and velocity
    newX[newIndex] = new double[3];
    newV[newIndex] = new double[3];

    // Initialize mass and position/velocity with the body's current state
    newM[newIndex] = mass[i];
    for (int dim = 0; dim < 3; dim++)
    {
      newX[newIndex][dim] = x[i][dim];
      newV[newIndex][dim] = v[i][dim];
    }

    // Merge other bodies into this one
    for (int j = 0; j < NumberOfBodies; j++)
    {
      if (merge[j] == i)
      { // Body j merges into body i
        double totalMass = newM[newIndex] + mass[j];
        for (int dim = 0; dim < 3; dim++)
        {
          newX[newIndex][dim] = (newX[newIndex][dim] * newM[newIndex] + x[j][dim] * mass[j]) / totalMass;
          newV[newIndex][dim] = (newV[newIndex][dim] * newM[newIndex] + v[j][dim] * mass[j]) / totalMass;
        }
        newM[newIndex] = totalMass;
      }
    }

    newIndex++;
  }
  // Free old arrays
  for (int i = 0; i < NumberOfBodies; i++)
  {
    delete[] x[i];
    delete[] v[i];
  }
  delete[] x;
  delete[] v;
  delete[] mass;

  // Update pointers to new arrays
  x = newX;
  v = newV;
  mass = newM;
  NumberOfBodies = newNBodies;
  delete[] merge;
}

void NBodySimulation::updateBody()
{

  timeStepCounter++;
  maxV = 0.0;
  minDx = std::numeric_limits<double>::max();

  // force1 = force along x direction
  // force2 = force along y direction
  // force3 = force along z direction

  for (int j = 0; j < NumberOfBodies; j++)
  {
    force1[j] = 0.0;
    force2[j] = 0.0;
    force3[j] = 0.0;
  }
  for (int j = 0; j < NumberOfBodies; j++)
  {

    for (int i = 0; i < NumberOfBodies; i++)
    {
      if (i == j)
        continue;

      // x,y,z forces acting on particle j by i.
      force1[j] += force_calculation(i, j, 0);
      force2[j] += force_calculation(i, j, 1);
      force3[j] += force_calculation(i, j, 2);
    }
  }
  // seperate blocks so that we are not changing positions of bodies
  // while still calculating interactions with other bodies

  for (int j = 0; j < NumberOfBodies; j++)
  {

    x[j][0] = x[j][0] + timeStepSize * v[j][0];
    x[j][1] = x[j][1] + timeStepSize * v[j][1];
    x[j][2] = x[j][2] + timeStepSize * v[j][2];

    v[j][0] = v[j][0] + timeStepSize * force1[j] / mass[j];
    v[j][1] = v[j][1] + timeStepSize * force2[j] / mass[j];
    v[j][2] = v[j][2] + timeStepSize * force3[j] / mass[j];
    double V = std::sqrt(v[j][0] * v[j][0] + v[j][1] * v[j][1] + v[j][2] * v[j][2]);
    if (maxV < V)
    {
      maxV = V;
    }
  }
  // check_collision();

  t += timeStepSize;
}

/**
 * Check if simulation has been completed.
 */
bool NBodySimulation::hasReachedEnd()
{
  return t > tFinal;
}

void NBodySimulation::takeSnapshot()
{
  if (t >= tPlot)
  {
    printParaviewSnapshot();
    printSnapshotSummary();
    tPlot += tPlotDelta;
  }
}

void NBodySimulation::openParaviewVideoFile()
{
  videoFile.open("paraview-output/result.pvd");
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\""
               " version=\"0.1\""
               " byte_order=\"LittleEndian\""
               " compressor=\"vtkZLibDataCompressor\">"
            << std::endl
            << "<Collection>";
}

void NBodySimulation::closeParaviewVideoFile()
{
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
  videoFile.close();
}

void NBodySimulation::printParaviewSnapshot()
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
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
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

void NBodySimulation::printSnapshotSummary()
{
  std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t=" << t
            << ",\t dt=" << timeStepSize
            << ",\t v_max=" << maxV
            << ",\t dx_min=" << minDx
            << std::endl;
}

void NBodySimulation::printSummary()
{
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
}
