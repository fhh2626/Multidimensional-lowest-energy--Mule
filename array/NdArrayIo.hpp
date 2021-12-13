#ifndef NDARRAYIO_HPP
#define NDARRAYIO_HPP

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "NdArray.hpp"
#include "pystring.h"
#include "pystring.cpp"

// usage:
//   // read a dat file into an NdArray
//   // the second parameter determines the data type
//   auto a = NdArray::readDat(file, 0.0);
//   // write and NdArray into a external file
//   NdArray::writeDat("file.txt", arr);
//
// note:
//   if this file is included, one must also include pystring
//

namespace NdArray {

    // read a datfile into a 2d NdArray
    template<typename T>
    NdArray<T> readDat(const std::string& file, T dummyVar) {

        static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value, "T must be a number");

        std::ifstream readFile;
        readFile.open(file, std::ios::in);
        if (!readFile.is_open()) {
            std::cerr << "file cannot open!" << std::endl;
            exit(1);
        }

        // all lines in memory
        std::vector<std::vector<std::string> > allLines;

        // used for parsing
        std::string line;
        std::vector<std::string> splitedLine;

        // read the file into memory
        while (getline(readFile, line)) {
            if (pystring::startswith(line,"#")) {
                continue;
            }
            pystring::split(line, splitedLine);
            if (splitedLine.size() == 0) {
                continue;
            }
            allLines.push_back(splitedLine);
        }

        // the shape of the matrix
        // the col num is determined by the first row
        std::vector<int> shape = {int(allLines.size()), int(allLines[0].size())};
        NdArray<T> data(shape);

        for (int i = 0; i < shape[0]; i++) {
            for (int j = 0; j < shape[1]; j++) {
                if (j < allLines[i].size()) {
                    data[{i, j}] = std::stod(allLines[i][j]);
                }
                else {
                    data[{i, j}] = 0;
                }
            }
        }

        readFile.close();

        return data;
    }

    // write an 2d NdArray to a datFile
    template<typename T>
    void writeDat(const std::string& file, NdArray<T> arr) {

        std::ofstream writeFile;
        writeFile.open(file, std::ios::out);
        if (!writeFile.is_open()) {
            std::cerr << "file cannot open!" << std::endl;
            exit(1);
        }

        auto shape = arr.getShape();
        for (int i = 0; i < shape[0]; i++) {
            for (int j = 0; j < shape[1]; j++) {
                writeFile << arr[{i,j}] << " ";
            }
            writeFile << "\n";
        }

        writeFile.close();
    }
}

#endif // NDARRAYIO_HPP
