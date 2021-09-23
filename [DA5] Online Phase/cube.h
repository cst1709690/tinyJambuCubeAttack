#include "components.h"
#ifndef CUBE_H
#define CUBE_H
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <random>
#include <cstdio>

class Cube
{
    public:
        /** constructor and destructor **/
        Cube();
        ~Cube();

        /** auxiliary functions **/
        bool readCubesFromFile();
        bool readCoefficientsFromFile();
        bool checkValues();

        /** main functions **/
        void startRHSGeneration();
        void generateRHSValue();

        /** printing functions **/
        void printReadCube();
        void printValueTable();

        void printCubeToFile(std::ofstream&);
        void printSuperpoliesToFile(std::ofstream&);
        void printRHSValueToFile(std::ofstream&);

    protected:

    private:
        /** Parameters for cube attack **/
        int cubesNo;
        int ltNo;
        int lc;
        int mlen;
        int adlen;
        int NROUND1;
        int NROUND2;
        int totalCubes;
        int falseCube;

        /** Variables used in computation **/
        uint8_t randomKey[16];
        uint32_t randomState[4];
        std::vector<std::vector<int>> superpolies;
        std::vector<uint8_t> rhsValues;
        std::vector<int> readCube;
        std::vector<int> bitExtractionShift;

        /** Input files used **/
        std::ifstream fileParameters;
        std::ifstream fileCubePassed;
        std::ifstream fileCoefficients;

        /** Output files used **/
        std::ofstream fileFilteredRHS;
        std::ofstream fileFilteredLHS;
        std::ofstream fileFilteredCube;

        std::ofstream fileFalseCube;
        std::ofstream fileFalseLHS;
        std::ofstream fileFalseRHS;
};

/** constructor **/
Cube::Cube()
{
    /** Getting parameters from the file **/
    std::string tempStr;
    std::stringstream tempSS;

    fileParameters.open("Parameters.txt");
    getline(fileParameters, tempStr);     // getting no. of cubes tested (not used)
    getline(fileParameters, tempStr);     // getting no. of tests (not used)

    getline(fileParameters, tempStr);       // getting cube size
    tempSS << tempStr;
    tempSS >> (this -> lc);
    tempSS.clear();

    getline(fileParameters, tempStr);       // getting plaintext size
    tempSS << tempStr;
    tempSS >> (this -> mlen);
    tempSS.clear();

    getline(fileParameters, tempStr);       // getting associated data size
    tempSS << tempStr;
    tempSS >> (this -> adlen);
    tempSS.clear();

    getline(fileParameters, tempStr);       // getting no. of rounds
    tempSS << tempStr;
    tempSS >> (this -> NROUND1);
    tempSS.clear();

    getline(fileParameters, tempStr);
    tempSS << tempStr;
    tempSS >> (this -> NROUND2);
    tempSS.clear();

    /** closing files used in constructors and opening other files for use **/
    fileParameters.close();

    /** Building the random key **/
    randomKey[0] = 0x00;
    randomKey[1] = 0x11;
    randomKey[2] = 0x22;
    randomKey[3] = 0x33;
    randomKey[4] = 0x44;
    randomKey[5] = 0x55;
    randomKey[6] = 0x66;
    randomKey[7] = 0x77;
    randomKey[8] = 0x88;
    randomKey[9] = 0x99;
    randomKey[10] = 0xAA;
    randomKey[11] = 0xBB;
    randomKey[12] = 0xCC;
    randomKey[13] = 0xDD;
    randomKey[14] = 0xEE;
    randomKey[15] = 0xFF;

    /** Building the random state **/
    randomState[0] = 0x33221100;
    randomState[1] = 0x77665544;
    randomState[2] = 0xBBAA9988;
    randomState[3] = 0xFFEEDDCC;

    /** cleaning up vectors and clearing constant **/
    readCube.clear();
    rhsValues.clear();
    superpolies.clear();
    bitExtractionShift.clear();

     /** counting the lines in successfulCube.txt **/
    int line = 0;
    fileCubePassed.open("successfulCubes.txt");
    while(getline(fileCubePassed, tempStr))
        line++;
    this -> totalCubes = line / 2;
    fileCubePassed.close();

    /** Opening input files **/
    fileCubePassed.open("successfulCubes.txt");
    fileCoefficients.open("LHS.txt");
    getline(fileCoefficients, tempStr);
    getline(fileCoefficients, tempStr);

    /** Opening output files **/
    fileFilteredRHS.open("filteredRHS.txt");
    fileFilteredCube.open("filteredCubes.txt");
    fileFilteredLHS.open("filteredSuperpolies.txt");

    fileFalseCube.open("falseCubes.txt");
    fileFalseLHS.open("falseSuperpolies.txt");
    fileFalseRHS.open("falseRHS.txt");
}

/** destructor **/
Cube::~Cube()
{
    free(randomKey);
}

/** auxiliary functions **/
bool Cube::readCubesFromFile()
{
    /** reading cube from file **/
    std::string fileCubeLine;
    std::string fileCoefficientLine;
    std::string cubeIndexStr;
    std::string coefficientIndexStr;

    /** clearing the vectors **/
    superpolies.clear();
    readCube.clear();
    rhsValues.clear();
    bitExtractionShift.clear();

    /** initialising the rhsValue vector **/
    for(int i = 0; i < 32; i++)
        rhsValues.push_back(-1);

    /** getting a line from the Cubes file **/
    std::getline(fileCubePassed, fileCubeLine);

    if(!fileCubePassed.eof() && fileCubeLine.size())
    {
        std::istringstream cubeSS(fileCubeLine);
        getline(cubeSS, cubeIndexStr, ' ');     // getting the empty string from the first space
        while(std::getline(cubeSS, cubeIndexStr, ' '))
        {
            std::stringstream cubeStoi(cubeIndexStr);
            int cubeIndexInt = 0;
            cubeStoi >> cubeIndexInt;
            readCube.push_back(cubeIndexInt);
        }
    }


    /** getting 32 lines from the coefficients file (to determine this cube passes which index of the keystream) **/
    if(!fileCoefficients.eof() && fileCubeLine.size())
    {
        for(int i = 0; i < 32; i++)
        {
            std::vector<int> tempCoefficientInts;

            /** getting a line for each 32 bit of keystream **/
            std::getline(fileCoefficients, fileCoefficientLine);
            std::istringstream coefficientsSS(fileCoefficientLine);

            /** processing coefficient **/
            while(std::getline(coefficientsSS, coefficientIndexStr, ' '))
            {
                std::stringstream coefficientStoi(coefficientIndexStr);
                int tempCoefficientInt;
                coefficientStoi >> tempCoefficientInt;
                tempCoefficientInts.push_back(tempCoefficientInt);
            }

            if(tempCoefficientInts[0] != -1)
                bitExtractionShift.push_back(i);

            superpolies.push_back(tempCoefficientInts);
        }

        getline(fileCoefficients, fileCoefficientLine);
    }


    if(!fileCubePassed.eof())
        return true;
    else
        return false;
}


/** main functions **/
void Cube::startRHSGeneration()
{
    int counter = 0;    // progress counter

    while(readCubesFromFile())
    {
        if(readCube.size())
        {
            std::cout << std::endl;
            std::cout << std::endl << "Overall Progress\t: " << ++counter << " / " << totalCubes;
            printReadCube();
            std::cout << std::endl << "Cube Size\t\t: " << lc;
            std::cout << std::endl << "Plaintext Length\t: " << mlen;
            std::cout << std::endl << "Associated Data Length\t: " << adlen;
            std::cout << std::endl << "Rounds\t\t\t: " << NROUND2;
            generateRHSValue();
            bool correctCube = checkValues();
            std::cout << std::endl << "Correct Cube\t\t: " << correctCube ? "True" : "False";
            printValueTable();
            if(correctCube)
            {
                printCubeToFile(fileFilteredCube);
                printSuperpoliesToFile(fileFilteredLHS);
                printRHSValueToFile(fileFilteredRHS);
            }
            else
            {
                printCubeToFile(fileFalseCube);
                printSuperpoliesToFile(fileFalseLHS);
                printRHSValueToFile(fileFalseRHS);
            }
        }
    }
}

void Cube::generateRHSValue()
{
    uint8_t plaintext[mlen] = {0};
    uint8_t ad[adlen] = {0};
    uint8_t tempCubeCombination;
    uint32_t tempKeystream = 0;
    uint32_t tempState[4] = {0};

    uint64_t combinations = 1ULL << lc;
    uint64_t trackCombinations = -1;

    for(int i = 0; i < combinations; i++)
    {
        trackCombinations++;
        for(int j = 0; j < lc; j++)
        {
            tempCubeCombination = (trackCombinations >> j) & 0x01;
            if(tempCubeCombination)
                ad[readCube[j] >> 3] |= (1 << (readCube[j] & 7));
            else
                ad[readCube[j] >> 3] &= (~(1 << (readCube[j] & 7)));
        }

        tempState[0] = randomState[0];
        tempState[1] = randomState[1];
        tempState[2] = randomState[2];
        tempState[3] = randomState[3];

        TinyJambu tinyJambu(NROUND1, NROUND2);
        tinyJambu.process_ad(randomKey, ad, adlen, tempState);
        uint32_t ks = tinyJambu.process_plaintext(randomKey, plaintext, mlen, tempState);

        tempKeystream ^= ks;
    }

    for(int i = 0; i < bitExtractionShift.size(); i++)
        rhsValues[bitExtractionShift[i]] = (tempKeystream & (1 << bitExtractionShift[i])) >> bitExtractionShift[i];
}

bool Cube::checkValues()
{
    for(int i = 0; i < bitExtractionShift.size(); i++)
    {
        uint8_t rhsValue = superpolies[bitExtractionShift[i]][0];
        if(superpolies[bitExtractionShift[i]].size() > 1)
        {
            for(int j = 1; j < superpolies[bitExtractionShift[i]].size(); j++)
            {
                int idx = superpolies[bitExtractionShift[i]][j];
                if(idx < 128)
                    rhsValue ^= ((randomKey[idx >> 3] & (1 << (idx & 7))) >> (idx & 7));
                else
                {
                    idx -= 128;
                    rhsValue ^= ((randomState[idx >> 5] & (1 << (idx & 31))) >> (idx & 31));
                }
            }

            if(rhsValue != rhsValues[bitExtractionShift[i]])
                return false;
        }
        else if(superpolies[bitExtractionShift[i]].size() == 1)
        {
            if(rhsValue != rhsValues[bitExtractionShift[i]])
                return false;
        }
    }

    return true;
}

/** printing functions **/
void Cube::printReadCube()
{
    std::vector<int>::iterator itr;

    std::cout << std::endl << "Cubes Read from File:\t: ";
    for(itr = readCube.begin(); itr != readCube.end(); itr++)
        std::cout << *itr << " ";
}

void Cube::printValueTable()
{
    std::cout << std::endl << "\t\tSuperpolies\t\tRHS Value";
    for(int i = 0; i < superpolies.size(); i++)
    {
        std::cout << std::endl << "State[" << 64 + i << "]:\t" << superpolies[i][0];
        for(int j = 1; j < superpolies[i].size(); j++)
            std::cout << " " << superpolies[i][j];
        std::cout << "\t\t\t" << rhsValues[i];
    }
}

void Cube::printCubeToFile(std::ofstream& theFile)
{
    theFile << std::endl;
    for(int i = 0; i < readCube.size(); i++)
        theFile << " " << readCube[i];
}

void Cube::printSuperpoliesToFile(std::ofstream& theFile)
{
    theFile << std::endl;
    for(int i = 0; i < superpolies.size(); i++)
    {
        for(int j = 0; j < superpolies[i].size(); j++)
            theFile << " " << superpolies[i][j];
        theFile << ", ";
    }
}

void Cube::printRHSValueToFile(std::ofstream& theFile)
{
    theFile << std::endl;
    for(int i = 0; i < rhsValues.size(); i++)
        theFile << " " << rhsValues[i];
}
#endif // CUBE_H
