#include "components.h"
#ifndef CUBE_H
#define CUBE_H
#include <vector>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <random>

class Cube
{
    public:
        /** constructor and destructor **/
        Cube(int, int, int, int, uint32_t, uint32_t);
        ~Cube();

        /** auxiliary functions **/
        void seed();
        void generateCube();

        /** main functions **/
        void startAttack();
        void linearityTestZero();
        bool linearityTest();
        void coefficientGeneration();

        /** printing functions **/
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
        int mlen;

        /** Parameters for initialisation rounds **/
        uint32_t NROUND1;
        uint32_t NROUND2;

        /** Constant (plaintext) **/
        uint8_t *plaintext;

        /** Variables used in computation **/
//        uint8_t *keyMatrix[128];   // key matrix including the zero key
        uint32_t constant;
        std::vector<int> generatedCube;
        std::vector<std::vector<uint8_t>> coefficients;
        std::vector<int> bitExtractionShift;

        /** Files used **/
        std::ofstream fileParameters;
        std::ofstream fileCubePassed;
        std::ofstream fileCoefficients;
};

/** constructor **/
Cube::Cube(int cubesNo = 0, int ltNo = 0, int lc = 0, int mlen = 0, uint32_t NROUND1 = 384, uint32_t NROUND2 = 1024)
{
    /** Initialising members **/
    this -> cubesNo = cubesNo;
    this -> ltNo = ltNo;
    this -> lc = lc;
    this -> mlen = mlen;
    this -> NROUND1 = NROUND1;
    this -> NROUND2 = NROUND2;

    /** setting up the plaintext **/
    plaintext = (uint8_t *) malloc(sizeof(uint8_t) * mlen);
    for(int i = 0; i < mlen; i++)
        plaintext[i] = 0;

    /** cleaning up vectors **/
    generatedCube.clear();
    coefficients.clear();
    bitExtractionShift.clear();

    /** Generating the parameters file **/
    fileParameters.open("Parameters.txt");
    fileParameters  << cubesNo      << std::endl
                    << ltNo         << std::endl
                    << lc           << std::endl
                    << mlen         << std::endl
                    << NROUND1      << std::endl
                    << NROUND2;
    fileParameters.close();
}

/** destructor **/
Cube::~Cube()
{
    /** freeing the keyMatrix **/
    free(plaintext);

    /** closing the files **/
    fileParameters.close();
    fileCubePassed.close();
    fileCoefficients.close();
}

/** auxiliary functions **/
void Cube::seed()
{
    srand((unsigned) time(0));
}

void Cube::generateCube()
{
    generatedCube.clear();

    /** for testing specific cube indices **/
    generatedCube.push_back(64);
    generatedCube.push_back(65);
    generatedCube.push_back(66);
    generatedCube.push_back(67);
    generatedCube.push_back(68);
    generatedCube.push_back(69);
    generatedCube.push_back(70);
    generatedCube.push_back(71);
    generatedCube.push_back(72);
    generatedCube.push_back(73);
    generatedCube.push_back(74);
    generatedCube.push_back(75);
    generatedCube.push_back(76);
    generatedCube.push_back(77);
    generatedCube.push_back(78);
    generatedCube.push_back(79);
    generatedCube.push_back(80);
    generatedCube.push_back(81);
    generatedCube.push_back(82);
    generatedCube.push_back(83);
    generatedCube.push_back(84);
    generatedCube.push_back(85);
    generatedCube.push_back(86);
    generatedCube.push_back(87);
    generatedCube.push_back(88);
    generatedCube.push_back(89);
    generatedCube.push_back(90);
    generatedCube.push_back(91);
    generatedCube.push_back(92);
    generatedCube.push_back(93);
    generatedCube.push_back(94);
    generatedCube.push_back(95);
    lc = generatedCube.size();

    /** random cube generation **/
//    int temp;
//    bool found = 1;
//    std::vector<int>::iterator testitr;
//    if(!(generatedCube.empty()))
//        generatedCube.clear();
//    while(generatedCube.size() < lc)
//    {
//        temp = rand() % 96;
//        testitr = find(generatedCube.begin(), generatedCube.end(), temp);
//        if(testitr == generatedCube.end())
//        {
//            found = 0;
//            generatedCube.push_back(temp);
//        }
//        else
//            continue;
//    }
//    sort(generatedCube.begin(), generatedCube.end());

    /** clearing the vectors and constant for new cube **/
    bitExtractionShift.clear();
    coefficients.clear();
    constant = 0;

    /** setting up bitExtractionShift vector elements of 0 to 31 **/
    for(int i = 0; i < 32; i++)
        bitExtractionShift.push_back(i);
}

/** main functions **/
void Cube::startAttack()
{
    fileCubePassed.open("successfulCubes.txt");
    fileCoefficients.open("LHS.txt");

    for(int i = 0; i < cubesNo; i++)
    {
        std::cout << std::endl << "Overall Progress\t\t\t\t\t: " << i + 1 << " / " << cubesNo << " cubes";
        generateCube();
        printGeneratedCube();
        std::cout << std::endl << "Cube Size\t\t\t\t\t\t: " << lc;
        std::cout << std::endl << "Number of Initialisation Rounds\t\t\t\t: " << NROUND1 << ", " << NROUND2;
        std::cout << std::endl << "Plaintext Length (bytes)\t\t\t\t: " << mlen;
        std::cout << std::endl;

        /** testing cubes by bypassing Linearity Test **/
//            coefficientGeneration();
//            printCoefficients();

        /** normal program execution **/
        int pass = 0;
        linearityTestZero();
        for(int j = 0; j < ltNo; j++)
        {
            std::cout << "\rLinearity Test Progress for " << ltNo << " Rounds\t\t\t: " << j + 1 << " / " << ltNo;
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
            std::cout << "Passed All Linearity Tests\t\t\t\t: " << ltNo << " times" << std::endl;
            printSuccessfulCubeToFile();
            coefficientGeneration();
            printCoefficients();
            printCoefficientsToFile();
        }
        else
            std::cout << "Passed Linearity Test\t\t\t\t: " << pass << " times" << std::endl;

        std::cout << std::endl;
    }
}

void Cube::linearityTestZero()
{
    /** Declaring component variables **/
    uint8_t *zeroKey;
    uint8_t *nonce;
    uint32_t *zeroState;

    /** Declaring temporary variables **/
    uint8_t tempCombination;
    uint32_t tempKeystream = 0;

    /** Declaring combination variables **/
    uint64_t combinations = 1ULL << lc;
    uint64_t trackCombinations = -1;

    /** Constructing the zero key and state **/
    zeroKey = (uint8_t *)malloc(sizeof(uint8_t) * 16);
    zeroState = (uint32_t *)malloc(sizeof(uint32_t) * 4);

    for(int i = 0; i < 16; i++)
        zeroKey[i] = 0;
    for(int i = 0; i < 4; i++)
        zeroState[i] = 0;

    /** Constructing the nonce **/
    nonce = (uint8_t *)malloc(sizeof(uint8_t) * 12);
    for(int i = 0; i < 12; i++)
        nonce[i] = 0;

    for(int i = 0; i < combinations; i++)
    {
        trackCombinations++;
        for(int j = 0; j < lc; j++)
        {
            tempCombination = (trackCombinations >> j) & 0x01;
            if(tempCombination)
                nonce[generatedCube[j] / 8] |= (1 << (generatedCube[j] % 8));
            else
                nonce[generatedCube[j] / 8] &= (~(1 << (generatedCube[j] % 8)));
        }

        TinyJambu tinyJambu(NROUND1, NROUND2);
        tinyJambu.initialization(zeroKey, nonce, zeroState);
        uint32_t ks = tinyJambu.process_plaintext(zeroKey, plaintext, mlen, zeroState);
        tempKeystream ^= ks;

        for(int j = 0; j < 4; j++)
            zeroState[j] = 0;
    }

    constant = tempKeystream;

    free(zeroKey);
    free(zeroState);
    free(nonce);
}

bool Cube::linearityTest()
{
    /** Declaring component variables **/
    uint8_t *keys[3];
    uint8_t *nonce;
    uint32_t *states[3];

    /** Declaring temporary variables **/
    uint8_t tempCombination;
    uint32_t tempInitialState[4];
    uint32_t tempKeystream[3] = {0};
    uint32_t results;

    /** Declaring combination variables **/
    uint64_t combinations = 1ULL << lc;    // Can only keep track of cubes of size 64
    uint64_t trackCombinations = -1;

    /** Declaring randomness variables **/
    std::minstd_rand0 generator(std::chrono::system_clock::now().time_since_epoch().count());

    /** Constructing random keys and states **/
    for(int i = 0; i < 3; i++)
    {
        keys[i] = (uint8_t *)malloc(sizeof(uint8_t) * 16);
        states[i] = (uint32_t *)malloc(sizeof(uint32_t) * 4);
    }
    for(int i = 0; i < 16; i++)
    {
        *(keys[0] + i) = (uint8_t) rand();
        *(keys[1] + i) = (uint8_t) rand();
        *(keys[2] + i) = *(keys[0] + i) ^ *(keys[1] + i);
    }
    for(int i = 0; i < 4; i++)
    {
        *(states[0] + i) = (generator() << 1) | (rand() & 1);
        *(states[1] + i) = (generator() << 1) | (rand() & 1);
        *(states[2] + i) = *(states[0] + i) ^ *(states[1] + i);
    }

    /** Used for testing specific set of keys and states **/
//    uint8_t keys[3][16] = {{0xAC, 0x67, 0xE5, 0xBE, 0x4B, 0x56, 0x87, 0x11, 0x9E, 0xE3, 0x60, 0x5F, 0x25, 0xC7, 0x9A, 0x0C},
//                           {0xE7, 0xD0, 0xF8, 0xDA, 0xD3, 0xAF, 0x8F, 0xE7, 0x43, 0x4D, 0xB7, 0xCB, 0xD8, 0x1D, 0x6B, 0x10},
//                           {0x4B, 0xB7, 0x1D, 0x64, 0x98    , 0xF9, 0x08, 0xF6, 0xDD, 0xAE, 0xD7, 0x94, 0xFD, 0xDA, 0xF1, 0x1C}};
//    uint32_t states[3][4] = {{0xB6162A72, 0xA1C27C0C, 0x04464A62, 0xEA62354D},
//                            {0x0A5356E5, 0xFD505458, 0x77039374, 0xE2E4C381},
//                            {0xBC457C97, 0x5C922854, 0x7345D916, 0x0886F6CC}};

    /** Constructing nonce **/
    nonce = (uint8_t *)malloc(sizeof(uint8_t) * 12);
    for(int i = 0; i < 12; i++)
        nonce[i] = 0;

    for(int i = 0; i < combinations; i++)
    {
        trackCombinations++;
        for(int j = 0; j < lc; j++)
        {
            tempCombination = (trackCombinations >> j) & 0x01;

            if(tempCombination)
                    nonce[generatedCube[j] / 8] |= (1 << (generatedCube[j] % 8));
            else
                    nonce[generatedCube[j] / 8] &= (~(1 << (generatedCube[j] % 8)));
        }

        for(int j = 0; j < 3; j++)
        {
            /** Saving the random state (cuz will be modified in the function later) **/
            for(int k = 0; k < 4; k++)
                tempInitialState[k] = *(states[j] + k);

            /** Performing the function and accumulating the keystream **/
            TinyJambu tinyJambu(NROUND1, NROUND2);
            tinyJambu.initialization(keys[j], nonce, tempInitialState);
            uint32_t ks = tinyJambu.process_plaintext(keys[j], plaintext, mlen, tempInitialState);
            tempKeystream[j] ^= ks;
        }
    }

    /** Determining which bit of the 32-bit keystream the cube passes */
    results = tempKeystream[0] ^ tempKeystream[1] ^ tempKeystream[2] ^ constant;

    /** Determining the cube passes which bit of keystream **/
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
    free(states[0]);
    free(states[1]);
    free(states[2]);

    if(bitExtractionShift.size())
        return true;
    else
        return false;
}

void Cube::coefficientGeneration()
{
    /** Declaring component variables **/
    uint8_t* nonce;

    /** Declaring temporary variables **/
    uint8_t tempCubeCombination;
    uint32_t tempInitialState[4];
    const uint32_t keystreamNo = 256;
    uint32_t tempKeystream[keystreamNo] = {0};
    std::vector<uint8_t> stateCoefficientTemp;

    /** Declaring combination variables **/
    uint64_t combinations = 1ULL << lc;    // Can only keep track of cubes of size 64
    uint64_t trackCombinations;

    /** Constructing nonce **/
    nonce = (uint8_t *) malloc(sizeof(uint8_t) * 12);
    for(int i = 0; i < 12; i++)
        nonce[i] = 0;

    /** Generating the 32 constants and push it into the coefficients 2D vector **/
    for(int i = 0; i < bitExtractionShift.size(); i++)
    {
        /** Clearing the vector for the current keystream bit **/
        stateCoefficientTemp.clear();

        /** Pushing constant **/
        stateCoefficientTemp.push_back((constant & (0x1 << bitExtractionShift[i])) >> bitExtractionShift[i]);
        coefficients.push_back(stateCoefficientTemp);
    }

    /** Computing the coefficients for the linear terms **/
    for(int i = 0; i < keystreamNo; i++)
    {
        std::cout << "\rCoefficient (Linear Terms) Generation Progress\t\t: " << i + 1 << " / " << keystreamNo;

        /** Setting up the key or state **/
        uint8_t key[16] = {0};
        uint32_t state[4] = {0};
        if(i < 128)
            key[i / 8] |= (1 << (i % 8));
        else
            state[(i - 128) / 32] |= (1 << ((i - 128) % 32));

        /** Iterating through all combinations of cube indices **/
        trackCombinations = -1;
        for(int j = 0; j < combinations; j++)
        {
            /** Updating the nonce **/
            trackCombinations++;
            for(int k = 0; k < lc; k++)
            {
                tempCubeCombination = (trackCombinations >> k) & 0x01;
                if(tempCubeCombination)
                    nonce[generatedCube[k] / 8] |= (1 << (generatedCube[k] % 8));
                else
                    nonce[generatedCube[k] / 8] &= (~(1 << (generatedCube[k] % 8)));
            }

            /** Creating a backup for the current value of state to be used in the next combination of the current cube **/
            for(int k = 0; k < 4; k++)
                tempInitialState[k] = state[k];

            /** Generating and accumulating keystream **/
            TinyJambu tinyJambu(NROUND1, NROUND2);
            tinyJambu.initialization(key, nonce, tempInitialState);
            uint32_t ks = tinyJambu.process_plaintext(key, plaintext, mlen, tempInitialState);

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
void Cube::printGeneratedCube()
{
    std::vector<int>::iterator itr;

    std::cout << std::endl << "Generated Cubes\t\t\t\t\t\t: ";
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
                    std::cout << (int) coefficients[currentPassedIndex][j] << " ";
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

    for(int i = 0; i < 32; i++)
    {
        fileCoefficients << std::endl;
        if(currentPassedIndex < bitExtractionShift.size())
        {
            if(i == bitExtractionShift[currentPassedIndex])
            {
                for(int j = 0; j < coefficients[currentPassedIndex].size(); j++)
                    fileCoefficients << (int) coefficients[currentPassedIndex][j] << " ";
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
