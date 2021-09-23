#include "Cube.h"

using namespace std;

int main()
{

    /** Testing cubes **/
    /*** Parameters ***/
    int cubesNo = 1;      // Number of cubes to be tested
    int lt = 100;           // Number of rounds of linearity tests
    int lc = 4;             // Cube size

    int NROUND1 = 384;     // Number of initialisation rounds for the lfsr and nfsr
    int NROUND2 = 1024;     // Number of initialisation rounds for the acc

    Cube cube(cubesNo, lt, lc, NROUND1, NROUND2);
    cube.seed();
    cube.startAttack();

    return 0;
}
