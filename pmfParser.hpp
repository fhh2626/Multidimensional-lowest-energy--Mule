#ifndef PMFPARSER_HPP
#define PMFPARSER_HPP

#include <iomanip>
#include <cassert>
#include <vector>

#include "array/NdArray.hpp"
#include "array/NdArrayIo.hpp"
#include "commonTools.h"

// parsing pmf files
// usage:
//   // read NAMD formatted PMF file
//   auto a = pmf<double>("file.pmf")
//   // read plain PMF file
//   auto a = pmf<double>("file.pmf",{-20,0},{0.2,0.1},{20,3})
//   // write NAMD formmatted PMF file
//   a.writePmfFile("file2.pmf")
//   // get data
//   a[{-20,0}]
//   a.getPmfData()
//   a.getLowerboundary()
//   a.getUpperboundary()
//   a.getWidth()
//   a.getShape()
//   a.getDimension()
//   a.RCToInternal()
//   a.internalToRC()
//

namespace pmfParser {

    // pmf (T=double) or count (T=int) data
    template <typename T>
    class pmf {

    public:

        // initialize the pmf using NAMD formatted file
        pmf(const std::string& pmfFile) {

            static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value, "T must be a kind of number");

            std::ifstream readFile;
            readFile.open(pmfFile, std::ios::in);
            if (!readFile.is_open()) {
                std::cerr << "file cannot be opened!" << std::endl;
                exit(1);
            }

            // used for parsing
            std::string line;
            std::vector<std::string> splitedLine;

            getline(readFile, line);
            if (!pystring::startswith(line,"#")) {
                std::cerr << "This is not an NAMD PMF file!" << std::endl;
                exit(1);
            }

            // dimension
            pystring::split(line, splitedLine);
            this->dimension = std::stoi(splitedLine[1]);

            // lb, width and ub
            this->lowerboundary = std::vector<double>(this->dimension);
            this->upperboundary = std::vector<double>(this->dimension);
            this->width = std::vector<double>(this->dimension);
            this->shape = std::vector<int>(this->dimension);
            for (int i = 0; i < this->dimension; i++) {
                getline(readFile, line);
                pystring::split(line, splitedLine);
                this->shape[i] = std::stoi(splitedLine[3]);
                this->width[i] = std::stod(splitedLine[2]);
                this->lowerboundary[i] = std::stod(splitedLine[1]) + 0.5 * this->width[i];
                this->upperboundary[i] = this->lowerboundary[i] + this->width[i] * (this->shape[i] - 1);
            }

            // reading data
            this->data = new NdArray::NdArray<T>(this->shape);
            std::vector<double> RCPosition(this->dimension);
            while (getline(readFile, line)) {
                pystring::split(line, splitedLine);
                if (splitedLine.size() == 0) {
                    continue;
                }
                for (int i = 0; i < this->dimension; i++) {
                    RCPosition[i] = std::stod(splitedLine[i]);
                }
                (*(this->data))[this->RCToInternal(RCPosition)] = std::stod(splitedLine[this->dimension]);
            }

            readFile.close();
        }

        // initialize the pmf given lb, ub, width
        pmf(
            const std::string& pmfFile,
            const std::vector<double>& lowerboundary,
            const std::vector<double>& width,
            const std::vector<double>& upperboundary
            ) {

            static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value, "T must be a kind of number");

            assert(lowerboundary.size() == width.size());
            assert(lowerboundary.size() == upperboundary.size());

            this->lowerboundary = lowerboundary;
            this->upperboundary = upperboundary;
            this->width = width;
            this->dimension = lowerboundary.size();

            // the shape of internal data
            this->shape = std::vector<int>(this->dimension, 0);
            for (int i = 0; i < this->dimension; i++) {
                // +1 means the boundaries are included
                this->shape[i] = int((upperboundary[i] - lowerboundary[i] + commonTools::accuracy) / width[i]) + 1;
            }
            this->data = new NdArray::NdArray<T>(this->shape);

            // read pmfFile into memory
            auto datFormatData = NdArray::readDat(pmfFile, T(0));
            auto datShape = datFormatData.getShape();
            // convert it into internal data
            std::vector<double> RCPosition;
            for (int row = 0; row < datShape[0]; row++) {
                RCPosition = std::vector<double>(this->dimension, 0);
                for (int col = 0; col < datShape[1]; col++) {
                    if (col < this->dimension) {
                        RCPosition[col] = datFormatData[{row, col}];
                    }
                    else {
                        (*(this->data))[this->RCToInternal(RCPosition)] = T(datFormatData[{row, col}]);
                        break;
                    }
                }
            }
        }

        // write internal data to a pmf file
        // in NAMD pmf format!
        // note: PBCs are not recorded! So they are zeroes!
        void writePmfFile(const std::string& file) const {

            std::ofstream writeFile;
            writeFile.open(file, std::ios::out);
            if (!writeFile.is_open()) {
                std::cerr << "file cannot be opened!" << std::endl;
                exit(1);
            }

            // the head of NAMD pmf file
            writeFile << std::setw(2) << "# " << this->dimension << "\n";
            for (int i = 0 ; i < this->dimension; i++) {
                writeFile << "# "
                          << std::setw(10) << this->lowerboundary[i] - 0.5 * this->width[i]
                          << std::setw(10) << this->width[i]
                          << std::setw(10) << this->shape[i] << " "
                          << 0 << "\n";
            }
            writeFile << "\n";

            // iterate over any dimension
            int n = 0;
            std::vector<int> loopFlag(this->dimension, 0);
            while (n >= 0) {

                auto RC = this->internalToRC(loopFlag);
                for (auto coor:RC) {
                    writeFile << commonTools::round(coor, commonTools::decimal_acc) << " ";
                }
                writeFile << (*(this->data))[loopFlag] <<"\n";

                // mimic an nD for loop
                n = this->dimension - 1;
                while (n >= 0) {
                    loopFlag[n] += 1;
                    if (loopFlag[n] > this->shape[n] - 1) {
                        loopFlag[n] = 0;
                        n--;
                        writeFile << "\n";
                    }
                    else {
                        break;
                    }
                }
            }

            writeFile.close();
        }

        // get the data (ndarray) of the pmf
        const NdArray::NdArray<T>& getPmfData() const {
            return *(this->data);
        }

        // get lowerboundary, upperboundary, width, shape and dimension
        const std::vector<double>& getLowerboundary() const {
            return this->lowerboundary;
        }

        const std::vector<double>& getUpperboundary() const {
            return this->upperboundary;
        }

        const std::vector<double>& getWidth() const {
            return this->width;
        }

        const std::vector<int>& getShape() const {
            return this->shape;
        }

        int getDimension() const {
            return this->dimension;
        }

        // operator[] to get the desired item at a given RCPosition
        // one may note that the var type of operator[] is extremely important
        const T& operator[] (const std::vector<double>& RCPosition) const {
            return (*(this->data))[this->RCToInternal(RCPosition)];
        }

        // get the desired item using internal RC
        // this will make parsing PMF easier
        const T& operator[] (const std::vector<int>& internalPosition) const {
            return (*(this->data))[internalPosition];
        }

        // convert external/real reaction coordinate into the internal coordinate
        std::vector<int> RCToInternal(const std::vector<double>& RCPosition) const {
            assert(RCPosition.size() == this->dimension);

            std::vector<int> internalPosition(this->dimension);
            for (int i = 0; i < this->dimension; i++) {
                internalPosition[i] = int((RCPosition[i] - this->lowerboundary[i] + commonTools::accuracy) / this->width[i]);
            }
            return internalPosition;
        }

        // convert internal coordinate into external/real reaction coordinate
        std::vector<double> internalToRC(const std::vector<int>& internalPosition) const {
            assert(internalPosition.size() == this->dimension);

            std::vector<double> RCPosition(this->dimension);
            for (int i = 0; i < this->dimension; i++) {
                RCPosition[i] = double((internalPosition[i] * this->width[i]) + this->lowerboundary[i]);
            }
            return RCPosition;
        }

        ~pmf() {
            delete this->data;
        }

    private:

        // the data (free energy) of the pmf
        NdArray::NdArray<T>* data;
        // lowerboundary, upperboundary, width, and dimension
        std::vector<double> lowerboundary;
        std::vector<double> upperboundary;
        std::vector<double> width;
        // the shape of internal data
        std::vector<int> shape;
        int dimension;
    };
}

#endif
