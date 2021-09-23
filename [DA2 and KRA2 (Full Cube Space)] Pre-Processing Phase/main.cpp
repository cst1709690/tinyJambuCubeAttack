#include "Cube.h"

using namespace std;

int main()
{
    /** Testing cubes **/
    /*** Parameters ***/
    int cubesNo = 500;      // Number of cubes to be tested
    int lt = 50;           // Number of rounds of linearity tests
    int lc = 15;             // Cube size

    int NROUND1 = 384;     // Number of initialisation rounds for the lfsr and nfsr
    int NROUND2 = 416;     // Number of initialisation rounds for the acc
    int NROUND3 = 416;

    Cube cube(cubesNo, lt, lc, NROUND1, NROUND2, NROUND3);
    cube.seed();
    cube.startAttack();

    return 0;
}
