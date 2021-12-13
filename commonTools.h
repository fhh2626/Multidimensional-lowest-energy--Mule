#ifndef COMMONTOOLS_H
#define COMMONTOOLS_H

#include <cmath>
#include <vector>

namespace commonTools {

    // the accuracy of float num calculation
    const double accuracy = 1e-8;
    const double decimal_acc = (-(log10(accuracy))) - 1;

    inline double round(double number, int digit) {
        return (std::round(number * pow(10, digit)) / pow(10, digit));
    }

    bool vectorEqual(const std::vector<double>& vec1, const std::vector<double>& vec2) {
        if (vec1.size() != vec2.size()) {
            return false;
        }
        for (int i = 0; i < vec1.size(); i++) {
            if (abs(vec1[i] - vec2[i]) > accuracy) {
                return false;
            }
        }
        return true;
    }

    bool vectorInVectorOfVector(const std::vector<int>& vec, const std::vector<std::vector<int> >& vecList) {
        for (auto& vecItem:vecList) {
            if (vec == vecItem) {
                return true;
            }
        }
        return false;
    }

    bool vectorInVectorOfVector(const std::vector<double>& vec, const std::vector<std::vector<double> >& vecList) {
        for (auto& vecItem:vecList) {
            if (vectorEqual(vec, vecItem)) {
                return true;
            }
        }
        return false;
    }


}

#endif // header guard

