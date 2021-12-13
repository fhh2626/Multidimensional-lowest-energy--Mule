// MUltidimensional Least Energy finder (MULE) v0.20 beta
// by Haohao Fu (fhh2626_at_gmail.com)
//
// Usage:
//    mule.exe config.ini
//
// In config.ini:
//    [mule]
//    directory             =   ./ref.pmf
//    lowerboundary         =   -20, 0              //(can be omitted if NAMD pmf is provided)
//    upperboundary         =    20, 3              //(can be omitted if NAMD pmf is provided)
//    width                 =   0.2, 0.1            //(can be omitted if NAMD pmf is provided)
//    initial               =   -20, 1.0
//    end                   =    20, 1.0
//    pbc                   =     0, 0
//    writeExploredPoints   =     0                 //(unnecessary, defalut=0)
//    target                =    20, 1.0, 0.1, 0.0  //(unnecessary, defines the targeted points and force constants)
//                                                  //(can define more than one targets)
//

#include <cstdlib>
#include <iostream>
#include <vector>

#include "pathFinder.hpp"
#include "pmfParser.hpp"
#include "array/pystring.h"
#include "ini/INIReader.h"

// read NAMD pmf file
const pmfParser::pmf<double>* readPMF(const std::string& fileName) {
    return new pmfParser::pmf<double>(fileName);
}

// read general pmf file
const pmfParser::pmf<double>* readPMF(
             const std::string& pmfFile,
             const std::vector<double>& lowerboundary,
             const std::vector<double>& width,
             const std::vector<double>& upperboundary
            ) {
    return new pmfParser::pmf<double>(pmfFile, lowerboundary, width, upperboundary);
}

// write data to a file
void writeData(const std::string& file, const std::vector<std::vector<double> >& points) {
    std::ofstream writeFile;
    writeFile.open(file, std::ios::out);
    if (!writeFile.is_open()) {
        std::cerr << "Cannot open " << file << std::endl;
        exit(1);
    }
    for (const auto& result:points) {
        for (const auto& item:result) {
            writeFile << item << " ";
        }
        writeFile << std::endl;
    }
    writeFile.close();
}

// write numbers to a file
void writeData(const std::string& file, const std::vector<double>& data) {
    std::ofstream writeFile;
    writeFile.open(file, std::ios::out);
    if (!writeFile.is_open()) {
        std::cerr << "Cannot open " << file << std::endl;
        exit(1);
    }
    for (const auto& item:data) {
        writeFile << item << std::endl;
    }
    writeFile.close();
}

// find optimized pathway
// return the total number of points explored
int findPathway(
                 const pmfParser::pmf<double>* pmfInfo,
                 const std::vector<double>& initialPoint,
                 const std::vector<double>& endPoint,
                 const std::vector<bool>& pbc,
                 const std::string& outputPrefix,
                 std::vector<std::vector<double> >& targetedPoints,
                 std::vector<std::vector<double> >& forceConstants,
                 bool writeExploredPoints = false
                 ) {
    auto pathFind = pathFinder::pathFinder(*pmfInfo, initialPoint, endPoint, pbc);
    std::vector<std::vector<double> > results;
    std::vector<double> energyResults;

    if (targetedPoints.size() != 0 && forceConstants.size() != 0) {
        pathFind.setTargetedPoints(targetedPoints, forceConstants);
        pathFind.Dijkstra(&pathFinder::pathFinder::manhattonPotential);
    }
    else {
        pathFind.Dijkstra();
    }

    pathFind.getResults(results, energyResults);

    std::string trajFile = outputPrefix + ".traj";
    std::string energyFile = outputPrefix + ".energy";

    // if one wants to write explored points
    std::vector<std::vector<double> > exploredPoints;
    std::string exploredPointsFile = outputPrefix + ".explored";
    if (writeExploredPoints) {
        pathFind.getExploredPoints(exploredPoints);
    }

    // write traj and energy
    writeData(trajFile, results);
    writeData(energyFile, energyResults);

    // write explored points
    if (writeExploredPoints) {
        writeData(exploredPointsFile, exploredPoints);
    }

    delete pmfInfo;
    pmfInfo = nullptr;

    return pathFind.getExploredPointNum();
}

// read input pars from config file
void readConfig(
                const std::string& file,
                std::string& pmfPath,
                std::vector<double>& lowerboundary,
                std::vector<double>& width,
                std::vector<double>& upperboundary,
                std::vector<double>& initialPoint,
                std::vector<double>& endPoint,
                std::vector<bool>& pbc,
                std::string& outputPrefix,
                std::vector<std::vector<double> >& targetedPoints,
                std::vector<std::vector<double> >& forceConstant,
                bool & writeExploredPoints
               ) {
    INIReader reader(file);
    if (reader.ParseError() != 0) {
        std::cerr << "Can't load ini file\n";
        exit(1);
    }

    pmfPath = reader.Get("mule", "directory", "");
    // these vars will be converted to std::vector
    auto tempLowerboundary = reader.Get("mule", "lowerboundary", "");
    auto tempUpperboundary = reader.Get("mule", "upperboundary", "");
    auto tempWidth = reader.Get("mule", "width", "");
    auto tempInitial = reader.Get("mule", "initial", "");
    auto tempEnd = reader.Get("mule", "end", "");
    auto tempPbc = reader.Get("mule", "pbc", "");
    auto tempTargetedPointsAndFC = reader.Get("mule", "target", "");

    writeExploredPoints = reader.GetBoolean("mule", "writeExploredPoints", false);

    std::vector<std::string> tempLowerboundaryStr, tempUpperboundaryStr, tempWidthStr;
    std::vector<std::string> tempInitialStr, tempEndStr, tempPbcStr, tempTargetedPointsAndFCStr;
    if (tempLowerboundary != "" && tempUpperboundary != "" && tempWidth != "") {
        pystring::split(tempLowerboundary, tempLowerboundaryStr, ",");
        pystring::split(tempUpperboundary, tempUpperboundaryStr, ",");
        pystring::split(tempWidth, tempWidthStr, ",");
        for (auto& item: tempLowerboundaryStr) lowerboundary.push_back(std::stod(item));
        for (auto& item: tempUpperboundaryStr) upperboundary.push_back(std::stod(item));
        for (auto& item: tempWidthStr) width.push_back(std::stod(item));
    }

    pystring::split(tempInitial, tempInitialStr, ",");
    pystring::split(tempEnd, tempEndStr, ",");
    pystring::split(tempPbc, tempPbcStr, ",");
    for (auto& item: tempInitialStr) initialPoint.push_back(std::stod(item));
    for (auto& item: tempEndStr) endPoint.push_back(std::stod(item));
    for (auto& item: tempPbcStr) pbc.push_back(std::stoi(item));

    // targeted points
    if (tempTargetedPointsAndFC != "") {
        std::vector<double> tempPoint, tempFC;
        pystring::split(tempTargetedPointsAndFC, tempTargetedPointsAndFCStr, ",");
        for (int i = 0; i < tempTargetedPointsAndFCStr.size() / initialPoint.size() / 2; i++) {
            tempPoint = {};
            tempFC = {};
            for (int j = 0; j < initialPoint.size(); j++) {
                tempPoint.push_back(std::stod(tempTargetedPointsAndFCStr[initialPoint.size() * 2 * i + j]));
                tempFC.push_back(std::stod(tempTargetedPointsAndFCStr[initialPoint.size() * 2 * i + initialPoint.size() + j]));
            }
            targetedPoints.push_back(tempPoint);
            forceConstant.push_back(tempFC);
        }
    }

    std::vector<std::string> tempOutputPrefix;
    pystring::rpartition(pmfPath, ".", tempOutputPrefix);
    outputPrefix = tempOutputPrefix[0];
}

int main(int argc, char* argv[]) {

    std::cout << "MUltidimensional Least Energy finder (MULE) v0.20 beta\n" << std::endl;

    if (argc < 2) {
        std::cerr << "Error, a config file must be provided!" << std::endl;
        exit(1);
    }
    std::string cfgFile(argv[1]);
    std::string pmfPath;
    std::vector<double> lowerboundary;
    std::vector<double> upperboundary;
    std::vector<double> width;
    std::vector<double> initialPoint;
    std::vector<double> endPoint;
    std::vector<bool> pbc;
    bool writeExploredPoints;
    std::string outputPrefix;
    std::vector<std::vector<double> > targetedPoints;
    std::vector<std::vector<double> > forceConstants;

    if (cfgFile.size() == 0) {
        std::cerr << "Error, a config file must be provided!" << std::endl;
        exit(1);
    }

    readConfig(
               cfgFile,
               pmfPath,
               lowerboundary,
               width,
               upperboundary,
               initialPoint,
               endPoint,
               pbc,
               outputPrefix,
               targetedPoints,
               forceConstants,
               writeExploredPoints
               );

    bool NAMDpmf = (lowerboundary.size() == 0 || upperboundary.size() == 0 || width.size() == 0);
    if (NAMDpmf) {
        std::cout << "Reading NAMD PMF file " << pmfPath << std::endl;
        std::cout << "Lowerboundary, upperboundary and width will be read from the PMF file!" << std::endl;
    }
    else {
        std::cout << "Reading plain PMF file " << pmfPath << std::endl;

        std::cout << "lowerboundary: " ;
        for (const auto& item: lowerboundary) std::cout << item << " ";
        std::cout << std::endl;

        std::cout << "upperboundary: ";
        for (const auto& item: upperboundary) std::cout << item << " ";
        std::cout << std::endl;

        std::cout << "width: ";
        for (const auto& item: width) std::cout << item << " ";
        std::cout << std::endl;
    }

    std::cout << "initial point: ";
    for (const auto& item: initialPoint) std::cout << item << " ";
    std::cout << std::endl;

    std::cout << "end point: ";
    for (const auto& item: endPoint) std::cout << item << " ";
    std::cout << std::endl;

    if (targetedPoints.size() != 0 && forceConstants.size() != 0) {
        std::cout << "Target points: " << std::endl;
        for (const auto& point: targetedPoints) {
            for (const auto& item: point) {
                std::cout << item << " ";
            }
            std::cout << std::endl;
        }
    }

    auto pmfInfo = NAMDpmf ? readPMF(pmfPath) : readPMF(pmfPath, lowerboundary, width, upperboundary);

    int exploredPointNum = findPathway(
                                       pmfInfo,
                                       initialPoint,
                                       endPoint,
                                       pbc,
                                       outputPrefix,
                                       targetedPoints,
                                       forceConstants,
                                       writeExploredPoints
                                      );

    std::cout << "Finished! See " << outputPrefix + ".traj" << " and " << outputPrefix + ".energy" << " for the results\n";
    std::cout << "A total of " << exploredPointNum << " points have been explored!\n";
    return 0;
}
