#include <qlat/qlat.h>
// #include <qlat/grid.h>

using namespace qlat;
using namespace std;

std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> mpi_coor({1, 1, 1, 8});

int main(int argc, char* argv[])
{
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h

  end();

  return 0;
}
