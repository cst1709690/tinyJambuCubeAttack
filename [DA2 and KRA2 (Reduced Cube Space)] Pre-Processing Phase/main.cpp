#include "Cube.h"

using namespace std;

int main()
{
    /*** Parameters ***/
    int cubesNo = 500;      // Does not matter since it is exhaustive cube search
    int lt = 100;           // Number of rounds of linearity tests
    int lc = 4;             // Cube size (less than 32)

    int NROUND1 = 384;      // Number of initialisation rounds for the key setup
    int NROUND2 = 1024;     // Number of initialisation rounds for the nonce setup
    int NROUND3 = 416;      // Number of initialisation rounds for the reduced encryption rounds

    Cube cube(cubesNo, lt, lc, NROUND1, NROUND2, NROUND3);
    cube.seed();
    cube.startAttack();

    return 0;
}
