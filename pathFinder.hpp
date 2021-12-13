#ifndef PATHFINDER_HPP
#define PATHFINDER_HPP

#include <cstdlib>
#include <cassert>
#include <vector>
#include <map>

#include "pmfParser.hpp"

// find the optimal pathway connecting two points on a pmf
// Usage:
//   // run calculation
//   auto path = pathFinder(pmfData, initialPoint, endPoint, pbc)
//   path.Dijkstra()
//   // one may want to add external manhatton potential
//   path.setTargetedPoints(
//                          {{19.5,2.2},{20.0,2.5}},
//                          {{1.0,1.0},{1.0,1.0}}
//                         )
//   path.Dijkstra(pathFinder::manhattonPotential)
//   // get results
//   std::vector<std::vector<double> > trajectory
//   std::vector<double> energyResults
//   path.getResults(trajectory, energyResults)
//   // get all the points explored during the process
//   std::vector<std::vector<double> > pointList
//   path.getExploredPoints(pointList)
//   // get the number of points explored
//   auto num = path.getExploredPointNum()
//

namespace pathFinder {

    // find the optimal pathway connecting two points on a PMF
    class pathFinder {

    public:

        // constructor
        pathFinder(
                   const pmfParser::pmf<double>& pmfData,
                   const std::vector<double>& initialPoint,
                   const std::vector<double>& endPoint,
                   const std::vector<bool>& pbc
                   ) {

            assert(initialPoint.size() == endPoint.size());
            assert(initialPoint.size() == pbc.size());
            assert(initialPoint.size() == pmfData.getDimension());

            this->pmfData = &pmfData;
            // internally, all the analyses are performed in the internal RC space
            this->lowerboundary = std::vector<int>(this->pmfData->getDimension(), 0);
            this->upperboundary = this->pmfData->getShape();
            // upperboundary is in fact shape - 1
            for (auto& b:this->upperboundary) {
                b -= 1;
            }
            this->width = std::vector<int>(this->pmfData->getDimension(), 1);
            this->initialPoint = this->pmfData->RCToInternal(initialPoint);
            this->endPoint = this->pmfData->RCToInternal(endPoint);
            this->pbc = pbc;
            this->dimension = this->pmfData->getDimension();

            this->openList.push_back(this->initialPoint);
        }

        // set targeted points and force constants
        // used in the A-star alg
        void setTargetedPoints (std::vector<std::vector<double> >& points, std::vector<std::vector<double> >& forceConst) {
            assert(points.size() == forceConst.size());
            for (int i = 0; i < points.size(); i++) {
                this->targetedPoints.push_back(this->pmfData->RCToInternal(points[i]));
                this->forceConstants.push_back(forceConst[i]);
            }
        }

        // run the Dijkstra alg
        void Dijkstra(double (pathFinder::*func)(const std::vector<int>& point) const = &pathFinder::defaultFunc) {

            // just a simple translation of the classical dijkstera alg
            while (this->openList.size() != 0) {
                auto p = this->popMin(this->openList, func);
                this->closeList.push_back(p);

                // the end point is found
                if (p == this->endPoint) {
                    break;
                }

                for (const auto& q:this->findAdjacentPoints(p)) {
                    if (commonTools::vectorInVectorOfVector(q, this->openList)) {
                        ; // pass
                    }
                    else if (!commonTools::vectorInVectorOfVector(q, this->closeList)) {
                        this->openList.push_back(q);
                        this->fatherPoint[q] = p;
                    }
                }
            }
        }

        // return the explored points of dijkstera calculation
        void getExploredPoints(std::vector<std::vector<double> > & pointList) const {

            if (this->closeList.size() == 0) {
                std::cerr << "Error, no information about results!\n";
                exit(1);
            }

            pointList = {};
            for (const auto &p:this->closeList) {
                pointList.push_back(this->pmfData->internalToRC(p));
            }
        }

        // how many points have been explored during the calculation
        int getExploredPointNum() const {
            if (this->closeList.size() == 0) {
                std::cerr << "Error, no information about results!\n";
                exit(1);
            }
            return this->closeList.size();
        }

        // return the minimum energy pathway of dijkstera calculation
        void getResults(std::vector<std::vector<double> >& trajectory, std::vector<double>& energyResults) {

            if (this->closeList.size() == 0) {
                std::cerr << "Error, no information about results!\n";
                exit(1);
            }

            std::vector<std::vector<int> > internalTrajectory = {};
            trajectory = {};
            energyResults = {};

            this->constructResults(this->endPoint, internalTrajectory);
            internalTrajectory.push_back(this->endPoint);

            for (const auto &p:internalTrajectory) {
                trajectory.push_back(this->pmfData->internalToRC(p));
            }

            for(const auto& p:internalTrajectory) {
                energyResults.push_back((*(this->pmfData))[p]);
            }
        }

        // below are h(x) functions used in the A-star alg
        inline double defaultFunc(const std::vector<int>& point) const {
            return 0;
        }

        // Manhatton potential
        double manhattonPotential(const std::vector<int>& point) const {
            if (this->targetedPoints.size() == 0) {
                return 0;
            }
            double energy = 0;
            int distance = 0;
            int d1, d2, d3;
            for (int i = 0; i < this->targetedPoints.size(); i++) {
                for (int j = 0; j < this->targetedPoints[i].size(); j++) {
                    distance = 0;
                    if (this->pbc[j] == false) {
                        distance = abs(point[j] - targetedPoints[i][j]);
                    }
                    else {
                        d1 = abs(point[j] - targetedPoints[i][j]);
                        d2 = (abs(point[j] - lowerboundary[j]) + abs(targetedPoints[i][j] - upperboundary[j]));
                        d3 = (abs(point[j] - upperboundary[j]) + abs(targetedPoints[i][j] - lowerboundary[j]));
                        distance = ((d1 < d2 ? d1 : d2) < d3) ? (d1 < d2 ? d1 : d2) : d3;
                    }
                    energy += distance * this->forceConstants[i][j];
                }
            }
            return energy;
        }

    private:

        // used in getResults
        std::vector<int> constructResults(std::vector<int> point, std::vector<std::vector<int> >& internalTrajectory) {
            if (this->fatherPoint.find(point) != this->fatherPoint.end()) {
                internalTrajectory.push_back(constructResults(this->fatherPoint[point], internalTrajectory));
            }
            return point;
        }

        // find the point that has the lowest energy in popList,
        // remove it from the list and return it
        std::vector<int> popMin(
                                std::vector<std::vector<int> >& popList,
                                double (pathFinder::*func)(const std::vector<int>& point) const = &pathFinder::defaultFunc
                               ) {

            assert(popList.size() != 0);

            // func(popList[0]) is h(x) in A-star alg
            // by default it is zero
            double minEnergy = (*pmfData)[popList[0]] + (this->*func)(popList[0]);
            int minPos = 0;

            for (int i = 0; i < popList.size(); i++) {
                if ((*pmfData)[popList[i]] + (this->*func)(popList[i]) < minEnergy) {
                    minEnergy = (*pmfData)[popList[i]] + (this->*func)(popList[i]);
                    minPos = i;
                }
            }

            auto minPoint = popList[minPos];
            popList.erase(popList.begin() + minPos);
            return minPoint;
        }

        // find the adjacent points of the input point,
        // return them as a list
        std::vector<std::vector<int> > findAdjacentPoints(const std::vector<int>& point) const {

            std::vector<std::vector<int> > adjacentPoints;

            for (int i = 0; i < this->dimension; i++) {
                // left side
                auto p = point;
                if (point[i] - this->width[i] >= this->lowerboundary[i] - commonTools::accuracy) {
                    p[i] = point[i] - this->width[i];
                }
                else {
                    if (this->pbc[i]) {
                        p[i] = this->upperboundary[i];
                    }
                    else {
                        p = {};
                    }
                }

                // right side
                auto q = point;
                if (point[i] + this->width[i] <= this->upperboundary[i] + commonTools::accuracy) {
                    q[i] = point[i] + this->width[i];
                }
                else {
                    if (this->pbc[i]) {
                        q[i] = this->lowerboundary[i];
                    }
                    else {
                        q = {};
                    }
                }

                if (p.size() != 0) {
                    adjacentPoints.push_back(p);
                }
                if (q.size() != 0) {
                    adjacentPoints.push_back(q);
                }
            }

            if (adjacentPoints.size() != 0) {
                return adjacentPoints;
            }
            else {
                std::cerr << "Error! No adjacent point is found!" << std::endl;
                exit(1);
            }
        }

        // the pmf data
        const pmfParser::pmf<double>* pmfData;
        std::vector<int> lowerboundary;
        std::vector<int> upperboundary;
        std::vector<int> width;
        // initial and end point
        std::vector<int> initialPoint;
        std::vector<int> endPoint;
        // whether periodic for each dimension
        std::vector<bool> pbc;
        // dimension of pmf
        int dimension;

        // below are vars that will be used in dijkstra/A* algs
        // record points' father point
        std::map<std::vector<int>, std::vector<int> > fatherPoint = {};
        // open and closeList in dijkstra/A* algs
        std::vector<std::vector<int> > openList = {};
        std::vector<std::vector<int> > closeList = {};

        // in the A* alg, one can define manhatton potential
        // based on targeted points and force constants
        std::vector<std::vector<int> > targetedPoints;
        std::vector<std::vector<double> > forceConstants;
    };
}

#endif // PATHFINDER_HPP
