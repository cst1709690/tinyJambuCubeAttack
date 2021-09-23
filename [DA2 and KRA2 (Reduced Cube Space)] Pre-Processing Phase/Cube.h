#include "components.h"
#ifndef CUBE_H
#define CUBE_H
#include <vector>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <sstream>

class Cube
{
    public:
        /** constructor and destructor **/
        Cube(int, int, int, int, int, int);
        ~Cube();

        /** auxiliary functions **/
        void seed();
        void initialiseCube();
        void readCubeSpace();
        void initialiseCubeSpace();

        /** main functions **/
        void startAttack();
        void buildExhaustiveCube(int, int, int);
        void preProcessingPhase();
        void linearityTestZero();
        bool linearityTest();
        void coefficientGeneration();

        /** printing functions **/
        void printCubeSpace();
        void printCubeSpaceToFile();
        void printGeneratedCube();
        void printSuccessfulCubeToFile();
        void printCoefficients();
        void printCoefficientsToFile();

    protected:

    private:
        /** Parameters for cube attack **/
        int cubesNo;
        int ltNo;
        int lc;
        int counter = 0;

        /** Parameters for the cube spaces **/
        int lcCubeSpace;
        int cubeSpaceNo;
        int testedCubeSpace = 0;
        std::vector<int> cubeSpace;

        /** Parameters for initialisation rounds **/
        int NROUND1;
        int NROUND2;
        int NROUND3;

        /** Variables used in computation **/
        uint32_t constant;
        uint8_t *keyMatrix[128];   // key matrix including the zero key
        std::vector<int> generatedCube;
        std::vector<std::vector<int>> coefficients;
        std::vector<int> bitExtractionShift;

        /** Files used **/
        std::ifstream fileCubeSpace;
        std::ofstream fileParameters;
        std::ofstream fileCubeSpacePassed;
        std::ofstream fileCubePassed;
        std::ofstream fileCoefficients;
};

/** constructor **/
Cube::Cube(int cubesNo = 0, int ltNo = 0, int lc = 0, int NROUND1 = 384, int NROUND2 = 1024, int NROUND3 = 1024)
{
    /** Initialising members **/
    this -> cubesNo = cubesNo;
    this -> ltNo = ltNo;
    this -> lc = lc;
    this -> NROUND1 = NROUND1;
    this -> NROUND2 = NROUND2;
    this -> NROUND3 = NROUND3;

    /** setting up the key matrix **/
    for(int i = 0; i < 128; i++)
    {
        keyMatrix[i] = (uint8_t *) malloc(sizeof(uint8_t) * 16);
        for(int j = 0; j < 16; j++)
        {
            if(i / 8 == j)
                *(keyMatrix[i] + j) = (uint8_t) 0x1 << (i % 8);
            else
                *(keyMatrix[i] + j) = 0x00;
        }
    }

    /** cleaning up vectors **/
    generatedCube.clear();
    coefficients.clear();
    bitExtractionShift.clear();

    /** Generating the parameters file **/
    fileParameters.open("Parameters.txt");
    fileParameters  << cubesNo      << std::endl
                    << ltNo         << std::endl
                    << lc           << std::endl
                    << NROUND1      << std::endl
                    << NROUND2      << std::endl
                    << NROUND3;
    fileParameters.close();
}

/** destructor **/
Cube::~Cube()
{
    /** freeing the keyMatrix **/
    for(int i = 0; i < 128; i++)
        free(keyMatrix[i]);

    /** closing the files **/
    fileCubeSpace.close();
    fileParameters.close();
    fileCubePassed.close();
    fileCubeSpacePassed.close();
    fileCoefficients.close();
}

/** auxiliary functions **/
void Cube::seed()
{
    srand((unsigned) time(0));
}

void Cube::initialiseCubeSpace()
{
    /** Set the number of cubes to test in the file **/
//    cubeSpaceNo = 1;

    /** Read number of cube spaces from file **/
    int line = 0;
    std::string tempStr;
    fileCubeSpace.open("cubeSpace.txt");
    while(getline(fileCubeSpace, tempStr))      // calculate total lines in the file
        line++;
    this -> cubeSpaceNo = line / 2;
    fileCubeSpace.close();

    /** Leave this here **/
    fileCubeSpace.open("cubeSpace.txt");        // open file to let the readCubeSpace function reads the cub
}

void Cube::readCubeSpace()
{
    cubeSpace.clear();

    /** Testing a cube space */
//    cubeSpace.push_back(64);
//    cubeSpace.push_back(65);
//    cubeSpace.push_back(66);
//    cubeSpace.push_back(67);
//    cubeSpace.push_back(68);
//    cubeSpace.push_back(69);
//    cubeSpace.push_back(70);
//    cubeSpace.push_back(71);
//    cubeSpace.push_back(72);
//    cubeSpace.push_back(73);
//    cubeSpace.push_back(74);
//    cubeSpace.push_back(75);
//    cubeSpace.push_back(76);
//    cubeSpace.push_back(77);
//    cubeSpace.push_back(78);
//    cubeSpace.push_back(79);
//    cubeSpace.push_back(80);
//    cubeSpace.push_back(81);
//    cubeSpace.push_back(82);
//    cubeSpace.push_back(83);
//    cubeSpace.push_back(84);
//    cubeSpace.push_back(85);
//    cubeSpace.push_back(86);
//    cubeSpace.push_back(87);
//    cubeSpace.push_back(88);
//    cubeSpace.push_back(89);
//    cubeSpace.push_back(90);
//    cubeSpace.push_back(91);
//    cubeSpace.push_back(92);
//    cubeSpace.push_back(93);
//    cubeSpace.push_back(94);
//    cubeSpace.push_back(95);
//    lcCubeSpace = cubeSpace.size();

    /** Read a cube space from file **/
    std::string tempStr;
    std::getline(fileCubeSpace, tempStr);
    tempStr = tempStr.substr(1, tempStr.size() - 1);

    /** Process the read cube **/
    std::istringstream cubeSS(tempStr);
    while(std::getline(cubeSS, tempStr, ' '))      // split the string of cube
    {
        /** Convert each cube index into integer **/
        std::stringstream cubeStoI(tempStr);
        int cubeIndexInt = 0;
        cubeStoI >> cubeIndexInt;
        cubeSpace.push_back(cubeIndexInt);
    }
    lcCubeSpace = cubeSpace.size();

    std::getline(fileCubeSpace, tempStr);           // get the blank line after the cube line
}

void Cube::initialiseCube()
{
    /** clearing the vectors and constant for new cube **/
    bitExtractionShift.clear();
    coefficients.clear();
    constant = 0;

    for(int i = 0; i < 32; i++)
        bitExtractionShift.push_back(i);
}

/** main functions **/
void Cube::startAttack()
{
    fileCubeSpacePassed.open("successfulCubeSpace.txt");
    fileCubePassed.open("successfulCubes.txt");
    fileCoefficients.open("LHS.txt");

    initialiseCubeSpace();
    testedCubeSpace = 0;
    while(testedCubeSpace < cubeSpaceNo)
    {
        counter = 0;
        readCubeSpace();

        int n = lcCubeSpace;
        int r = lc;
        buildExhaustiveCube(n, 0, r);

        testedCubeSpace++;
    }
}

void Cube::buildExhaustiveCube(int n, int left, int r)
{
    if(r == 0)
    {
        preProcessingPhase();
        return;
    }

    for(int i = left; i < n; i++)
    {
        generatedCube.push_back(cubeSpace[i]);
        buildExhaustiveCube(n, i + 1, r - 1);
        generatedCube.pop_back();
    }
}

void Cube::preProcessingPhase()
{
    std::cout << std::endl << "Cube Space Progress\t\t\t: " << testedCubeSpace + 1 << " / " << cubeSpaceNo << " cube spaces";
    printCubeSpace();
    std::cout << std::endl << "Cube Space Size\t\t\t\t: " << lcCubeSpace;
    std::cout << std::endl << "Testing Cube No.\t\t\t: " << ++counter;
    initialiseCube();
    printGeneratedCube();
    std::cout << std::endl << "Cube Size\t\t\t\t: " << lc;
    std::cout << std::endl << "Number of Initialisation Rounds\t\t: " << NROUND1 << ", " << NROUND2 << ", " << NROUND3;
    std::cout << std::endl;

    int pass = 0;
    linearityTestZero();
    for(int j = 0; j < ltNo; j++)
    {
        std::cout << "\rLinearity Test Progress for " << ltNo << " Rounds\t: " << j + 1 << " / " << ltNo;
        if(linearityTest())
        {
            pass++;
            continue;
        }
        else
            break;
    }
    std::cout << std::endl;
    if(pass == ltNo)
    {
        std::cout << "Passed All Linearity Test\t\t: " << ltNo << " times" << std::endl;
        printCubeSpaceToFile();
        printSuccessfulCubeToFile();
        coefficientGeneration();
        printCoefficients();
        printCoefficientsToFile();
    }
    else
        std::cout << "Passed Linearity Test\t\t\t: " << pass << " times" << std::endl;

    std::cout << std::endl;
}

void Cube::linearityTestZero()
{
    uint8_t *zeroKey;
    uint8_t *nonce;
    uint8_t temp;

    uint32_t tempKeystream = 0;

    uint64_t combinations = 1ULL << lc;
    uint64_t trackCombinations = -1;

    /** Constructing and initialising zeroKey **/
    zeroKey = (uint8_t *)malloc(sizeof(uint8_t) * 16);
    for(int i = 0; i < 16; i++)
        zeroKey[i] = 0;

    /** Constructing and initialising nonce **/
    nonce = (uint8_t *)malloc(sizeof(uint8_t) * 12);
    for(int i = 0; i < 12; i++)
        nonce[i] = 0;

    for(int i = 0; i < combinations; i++)
    {
        trackCombinations++;
        for(int j = 0; j < lc; j++)
        {
            temp = (trackCombinations >> j) & 0x01;
            if(temp)
                    nonce[generatedCube[j] / 8] |= (1 << (generatedCube[j] % 8));
            else
                    nonce[generatedCube[j] / 8] &= (~(1 << (generatedCube[j] % 8)));
        }

        TinyJambu tinyJambu(NROUND1, NROUND2, NROUND3);
        uint32_t ks = tinyJambu.initialization(zeroKey, nonce, tinyJambu.state);
        tempKeystream ^= ks;
    }

    constant = tempKeystream;

    free(zeroKey);
    free(nonce);
}

bool Cube::linearityTest()
{
    uint8_t *keys[3];
    uint8_t *nonce;
    uint8_t temp;
    uint32_t tempKeystream[3] = {0};
    uint32_t results;

    uint64_t combinations = 1ULL << lc;    // Can only keep track of cubes of size 64
    uint64_t trackCombinations = -1;

    /** Constructing and initialising the three keys **/
    for(int i = 0; i < 4; i++)
        keys[i] = (uint8_t *)malloc(sizeof(uint8_t) * 16);
    for(int i = 0; i < 16; i++)
    {
        *(keys[0] + i) = (uint8_t) rand();
        *(keys[1] + i) = (uint8_t) rand();
        *(keys[2] + i) = *(keys[0] + i) ^ *(keys[1] + i) ^ 0;
    }

    /** Constructing and initialising the nonce **/
    nonce = (uint8_t *)malloc(sizeof(uint8_t) * 12);
    for(int i = 0; i < 12; i++)
        nonce[i] = 0;

    for(int i = 0; i < combinations; i++)
    {
        trackCombinations++;
        for(int j = 0; j < lc; j++)
        {
            temp = (trackCombinations >> j) & 0x01;
            if(temp)
                    nonce[generatedCube[j] / 8] |= (1 << (generatedCube[j] % 8));
            else
                    nonce[generatedCube[j] / 8] &= (~(1 << (generatedCube[j] % 8)));
        }

        for(int j = 0; j < 3; j++)
        {
            TinyJambu tinyJambu(NROUND1, NROUND2, NROUND3);
            uint32_t ks = tinyJambu.initialization(keys[j], nonce, tinyJambu.state);
            tempKeystream[j] ^= ks;
        }
    }

    /** Looking at which bit of the 32-bit keystream the cube passes */
    results = tempKeystream[0] ^ tempKeystream[1] ^ tempKeystream[2] ^ constant;
    for(int i = 0; i < bitExtractionShift.size(); i++)
    {
        if(results & (1 << bitExtractionShift[i]))
        {
            std::vector<int>::iterator position = find(bitExtractionShift.begin(), bitExtractionShift.end(), bitExtractionShift[i]);
            if(position != bitExtractionShift.end())
            {
                bitExtractionShift.erase(position);
                i--;        // To compensate for the reduced size
            }
        }
    }

    free(nonce);
    free(keys[0]);
    free(keys[1]);
    free(keys[2]);

    if(bitExtractionShift.size())
        return true;
    else
        return false;
}

void Cube::coefficientGeneration()
{
    uint8_t *nonce;
    uint8_t temp;

    const uint32_t keystreamNo = 128;
    uint32_t tempKeystream[keystreamNo] = {0};

    uint64_t combinations = 1ULL << lc;    // Can only keep track of cubes of size 64
    uint64_t trackCombinations;
    std::vector<int> stateCoefficientTemp;

    /** initialise nonce **/
    nonce = (uint8_t *) malloc(sizeof(uint8_t) * 12);
    for(int i = 0; i < 12; i++)
        nonce[i] = 0x00;

    /** Generating the 32 constants and push it into the coefficients 2D vector **/
    for(int i = 0; i < bitExtractionShift.size(); i++)
    {
        /** Clearing the vector for the current keystream bit **/
        stateCoefficientTemp.clear();

        /** Pushing constant **/
        stateCoefficientTemp.push_back((constant & (0x1 << bitExtractionShift[i])) >> bitExtractionShift[i]);
        coefficients.push_back(stateCoefficientTemp);
    }

    for(int i = 0; i < keystreamNo; i++)
    {
        std::cout << "\rCoefficient Generation Progress\t\t: " << i + 1 << " / " << keystreamNo;
        trackCombinations = -1;
        for(int j = 0; j < combinations; j++)
        {
            trackCombinations++;
            for(int k = 0; k < lc; k++)
            {
                temp = (trackCombinations >> k) & 0x01;
                if(temp)
                    nonce[generatedCube[k] / 8] |= (1 << (generatedCube[k] % 8));
                else
                    nonce[generatedCube[k] / 8] &= (~(1 << (generatedCube[k] % 8)));
            }

            TinyJambu tinyJambu(NROUND1, NROUND2, NROUND3);
            uint32_t ks = tinyJambu.initialization(keyMatrix[i], nonce, tinyJambu.state);
            tempKeystream[i] ^= ks;
        }

        /** Iterate through and extract all the bits that pass the tests of the current accumulated keystream **/
        for(int j = 0; j < bitExtractionShift.size(); j++)
        {
            uint8_t currentKeystreamBit = (tempKeystream[i] & (0x1 << bitExtractionShift[j])) >> bitExtractionShift[j];
            if(coefficients[j][0] ^ currentKeystreamBit)
                coefficients[j].push_back(i);
        }
    }
    std::cout << std::endl;

    free(nonce);
}

/* printing functions */
void Cube::printCubeSpace()
{
    std::vector<int>::iterator itr;

    std::cout << std::endl << "Cube Space\t\t\t\t: ";
    for(itr = cubeSpace.begin(); itr != cubeSpace.end(); itr++)
        std::cout << *itr << " ";
}

void Cube::printCubeSpaceToFile()
{
    for(int i = 0; i < cubeSpace.size(); i++)
        fileCubeSpacePassed << " " << cubeSpace[i];
    fileCubeSpacePassed << std::endl << std::endl;
}

void Cube::printGeneratedCube()
{
    std::vector<int>::iterator itr;

    std::cout << std::endl << "Generated Cubes\t\t\t\t: ";
    for(itr = generatedCube.begin(); itr != generatedCube.end(); itr++)
        std::cout << *itr << " ";
}

void Cube::printSuccessfulCubeToFile()
{
    for(int i = 0; i < generatedCube.size(); i++)
        fileCubePassed << " " << generatedCube[i];
    fileCubePassed << std::endl << std::endl;
}

void Cube::printCoefficients()
{
    int currentPassedIndex = 0;
    int i = 0;

    for(int i = 0; i < 32; i++)
    {
        std::cout << std::endl << "Coefficients for state[" << 64 + i << "]: ";
        if(currentPassedIndex < bitExtractionShift.size())
        {
            if(i == bitExtractionShift[currentPassedIndex])
            {
                for(int j = 0; j < coefficients[currentPassedIndex].size(); j++)
                    std::cout << coefficients[currentPassedIndex][j] << " ";
                currentPassedIndex++;
            }
            else
                std::cout << -1;
        }
        else
            std::cout << -1;
    }
}

void Cube::printCoefficientsToFile()
{
    int currentPassedIndex = 0;

    fileCoefficients << std::endl;
    for(int i = 0; i < 32; i++)
    {
        fileCoefficients << std::endl;
        if(currentPassedIndex < bitExtractionShift.size())
        {
            if(i == bitExtractionShift[currentPassedIndex])
            {
                for(int j = 0; j < coefficients[currentPassedIndex].size(); j++)
                    fileCoefficients << coefficients[currentPassedIndex][j] << " ";
                currentPassedIndex++;
            }
            else
                fileCoefficients << -1;
        }
        else
            fileCoefficients << -1;
    }
}
#endif // CUBE_H
