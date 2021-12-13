#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include <cassert>
// #define NDEBUG

#include <iostream>
#include <vector>

// a lightweighted n-dimensional array library
// Haohao Fu (fhh2626@gmail.com)
// version 0.11 beta
//
// usage:
//   // initialize by shape (vector<int>) and default value (int, default 0)
//   NdArray<int> arr({5,4}, 1);
//   // copy constructor
//   NdArray<double> arr2 = arr;
//   // calculation
//   std::cout << arr2 * 5 + 1;
//   // reshape
//   // one must understand the mechanism of the reshape function
//   // internal data are stored in a 1d array
//   // NdArray.reshape() simply change the arrangement of the 1d array, for instance
//   // [1,2,3,4,5,6] (shape{6}) -> [[1,2,3],[4,5,6]] (shape{2,3}) -> [[1,2],[3,4],[5,6]] (shape{3,2})
//   arr.reshape({4,5});
//   // return the max/min value
//   arr.max()
//   // get other information of the array
//   arr.getShape()
//   arr.getTotalSize()
//   arr.getCArray()
//

namespace NdArray {

    // NdArray class
    template <typename T>
    class NdArray {
    public:
        // default constructor of NdArray
        NdArray(const std::vector<int>& shape, T defaultValue = 0) {

            static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value, "T must be a kind of number");

            // the shape of the tensor
            this->shape = shape;

            // calculate the length of the internal array
            this->totalSize = 1;
            for (auto i:shape) {
                this->totalSize *= i;
            }

            assert(this->totalSize != 0);

            // get memory
            this->data = new T[totalSize];
            // initialize the new array
            for (int i = 0; i < this->totalSize; i++) {
                this->data[i] = defaultValue;
            }
        }

        // copy constructor
        NdArray(const NdArray& arr) {

            this->totalSize = arr.getTotalSize();
            this->shape = arr.getShape();
            this->data = new T[this->totalSize];
            for (int i = 0; i < this->totalSize; i++) {
                this->data[i] = arr.getCArray()[i];
            }
        }

        // type change constructor
        template <typename U>
        NdArray(const NdArray<U>& arr) {

            this->totalSize = arr.getTotalSize();
            this->shape = arr.getShape();
            this->data = new T[this->totalSize];
            for (int i = 0; i < this->totalSize; i++) {
                this->data[i] = T(arr.getCArray()[i]);
            }
        }

        // operator[] to get the desired item at a given pos
        T& operator[] (const std::vector<int>& pos) {
            return data[position(pos)];
        }

        // operator +
        NdArray operator+ (const NdArray& arr) const {
            assert(this->shape == arr.shape);
            NdArray sumArr = arr;
            for (int i = 0; i < arr.totalSize; i++) {
                sumArr.data[i] += this->data[i];
            }
            return sumArr;
        }

        NdArray operator+ (const T& num) const {
            NdArray sumArr = *this;
            for (int i = 0; i < this->totalSize; i++) {
                sumArr.data[i] += num;
            }
            return sumArr;
        }

        friend NdArray operator+ (const T& num, const NdArray& arr) {
            return arr + num;
        }

        // operator +=
        NdArray operator+= (const NdArray& arr) {
            assert(this->shape == arr.shape);
            for (int i = 0; i < this->totalSize; i++) {
                this->data[i] += arr.data[i];
            }
            return *this;
        }

        NdArray& operator+= (const T& num) {
            for (int i = 0; i < this->totalSize; i++) {
                this->data[i] += num;
            }
            return *this;
        }

        // operator -
        NdArray operator- (const NdArray& arr) const {
            assert(this->shape == arr.shape);
            NdArray sumArr = *this;
            for (int i = 0; i < arr.totalSize; i++) {
                sumArr.data[i] -= arr.data[i];
            }
            return sumArr;
        }

        NdArray operator- (const T& num) const {
            NdArray sumArr = *this;
            for (int i = 0; i < this->totalSize; i++) {
                sumArr.data[i] -= num;
            }
            return sumArr;
        }

        // operator *
        NdArray operator* (const NdArray& arr) const {
            assert(this->shape == arr.shape);
            NdArray sumArr = *this;
            for (int i = 0; i < arr.totalSize; i++) {
                sumArr.data[i] *= arr.data[i];
            }
            return sumArr;
        }

        NdArray operator* (const T& num) const {
            NdArray sumArr = *this;
            for (int i = 0; i < this->totalSize; i++) {
                sumArr.data[i] *= num;
            }
            return sumArr;
        }

        friend NdArray operator* (const T& num, const NdArray& arr) {
            return arr * num;
        }

        // operator *=
        NdArray operator*= (const NdArray& arr) {
            assert(this->shape == arr.shape);
            for (int i = 0; i < this->totalSize; i++) {
                this->data[i] *= arr.data[i];
            }
            return *this;
        }

        NdArray& operator*= (const T& num) {
            for (int i = 0; i < this->totalSize; i++) {
                this->data[i] *= num;
            }
            return *this;
        }

        // operator /
        NdArray operator/ (const NdArray& arr) const {
            assert(this->shape == arr.shape);
            NdArray sumArr = *this;
            for (int i = 0; i < arr.totalSize; i++) {
                sumArr.data[i] /= arr.data[i];
            }
            return sumArr;
        }

        NdArray operator/ (const T& num) const {
            NdArray sumArr = *this;
            for (int i = 0; i < this->totalSize; i++) {
                sumArr.data[i] /= num;
            }
            return sumArr;
        }

        // operator %
        NdArray operator% (const NdArray& arr) const {
            assert(this->shape == arr.shape);
            NdArray sumArr = *this;
            for (int i = 0; i < arr.totalSize; i++) {
                sumArr.data[i] %= arr.data[i];
            }
            return sumArr;
        }

        NdArray operator% (const T& num) const {
            NdArray sumArr = *this;
            for (int i = 0; i < this->totalSize; i++) {
                sumArr.data[i] %= num;
            }
            return sumArr;
        }

        // operator <<
        friend std::ostream& operator<< (std::ostream& out, const NdArray& arr) {

            // iterate over any dimension
            int n = 0;
            std::vector<int> loopFlag(arr.shape.size(), 0);
            while (n >= 0) {

                out << arr.data[arr.position(loopFlag)] << " ";

                // mimic an nD for loop
                n = arr.shape.size() - 1;
                while (n >= 0) {
                    loopFlag[n] += 1;
                    if (loopFlag[n] > arr.shape[n] - 1) {
                        loopFlag[n] = 0;
                        n--;
                        out << "\n";
                    }
                    else {
                        break;
                    }
                }
            }
            return out;
        }

        // return the max value
        T maxValue() const {
            T maxN = this->data[0];
            for (int i = 0; i < this->totalSize; i++) {
                if (this->data[i] > maxN) {
                    maxN = data[i];
                }
            }
            return maxN;
        }

        // return the min value
        // return the max value
        T minValue() const {
            T minN = this->data[0];
            for (int i = 0; i < this->totalSize; i++) {
                if (this->data[i] < minN) {
                    minN = data[i];
                }
            }
            return minN;
        }

        // reshape the nd array
        void reshape(const std::vector<int>& newShape) {
            // calculate the length corresponding the new shape
            int totalSize = 1;
            for (auto i:newShape) {
                totalSize *= i;
            }

            assert(this->totalSize = totalSize);

            this->shape = newShape;
        }

        // return the totol size of the ndArray
        int getTotalSize() const {
            return this->totalSize;
        }

        // return the shape of the ndArray
        const std::vector<int>& getShape() const {
            return this->shape;
        }

        // get the C-style Array
        const T* const getCArray() const {
            return this->data;
        }

        // default destructor
        ~NdArray() {
            delete[] data;
        }

    private:

        // find the real position of an item based on a vector
        int position(const std::vector<int>& pos) const {
            assert(pos.size() == this->shape.size());

            int realPos = 0;
            int tempSum = 0;
            for (int i = 0; i < pos.size(); i++) {
                tempSum = pos[i];
                for (int j = i + 1; j < pos.size(); j++) {
                    tempSum *= this->shape[j];
                }
                realPos += tempSum;
            }
            return realPos;
        }

        // pointers to data
        T* data;
        // the total number of items
        int totalSize;
        // the shape of the NdArray
        std::vector<int> shape;
    };
}

#endif // NDARRAY_HPP


