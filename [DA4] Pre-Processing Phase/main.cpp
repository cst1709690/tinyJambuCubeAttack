#include "Cube.h"

using namespace std;

int main()
{
    /** Testing cubes **/
    /*** Parameters ***/
    int cubesNo = 1;      // Number of cubes to be tested
    int lt = 50;           // Number of rounds of linearity tests
    int lc = 32;             // Cube size
    int mlen = 4;           // Plaintext length in bytes

    int NROUND1 = 32;     // Number of initialisation rounds for the lfsr and nfsr
    int NROUND2 = 437;     // Number of initialisation rounds for the acc

    Cube cube(cubesNo, lt, lc, mlen, NROUND1, NROUND2);
    cube.seed();
    cube.startAttack();

    return 0;
}
