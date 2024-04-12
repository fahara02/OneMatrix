#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include <iomanip>
#include <limits>
#include <sstream>
// namespace MAT {

struct divide_by_zero : public std::exception {
  virtual const char* what() const throw() {
    return "division by zero occured!";
  }
};
struct bad_size : public std::exception {
  virtual const char* what() const throw() {
    return "matrix row or column size is incompatible";
  }
};
struct incompatible_size : public std::exception {
  virtual const char* what() const throw() {
    return "matrix/matricies not compatible sizes";
  }
};
struct bad_argument : public std::exception {
  virtual const char* what() const throw() {
    return "invalid argument to function";
  }
};

struct out_of_range : public std::exception {
  virtual const char* what() const throw() {
    return "row or column index is out of range";
  }
};
struct row_outbound : public out_of_range {
  const char* what() const throw() {
    return "indexing of matrix row is out of bound";
  }
};
struct col_outbound : public out_of_range {
  const char* what() const throw() {
    return "indexing of matrix column is out of bound";
  }
};
struct too_many_rows : public bad_argument {
  const char* what() const throw() { return "only one row at a time"; }
};

struct too_many_cols : public bad_argument {
  const char* what() const throw() { return "only one column at a time"; }
};

struct bad_rows : public bad_size {
  const char* what() const throw() {
    return "new row vector column size do not match ";
  }
};

struct bad_cols : public bad_size {
  const char* what() const throw() {
    return "new column vector row size do not match ";
  }
};
struct bad_row_match : public bad_size {
  const char* what() const throw() {
    return "new vector's row size dont match existing row ";
  }
};

struct bad_col_match : public bad_size {
  const char* what() const throw() {
    return "new vector's column size dont match existing column ";
  }
};
struct not_square : public bad_size {
  const char* what() const throw() { return "matrix must be square"; }
};
struct not_invertible : public std::exception {
  const char* what() const throw() { return "matrix is not invertible"; }
};
struct not_solvable : public std::exception {
  const char* what() const throw() { return "System is not solvable"; }
};

enum class VectorType { RowVector, ColumnVector };
enum class Orientation { Row, Column };

// template<typename T>
// struct Solution {
//     SolutionType type;
//     std::vector<T> values;
//     std::vector<std::string> linearComb;
// };

//enum class SliceType { Row, Column };

#define EPSILON 0.0000000001
#define EQUAL(a, b) (abs((a) - (b)) < EPSILON)

#define INDEX(r, c) ((m.mCols * (r)) + (c))

#define for_index(i, dim) for (int i = 0; i < (dim); ++i)

#define for_ij(r, c)               \
  for (size_t i = 0; i < (r); ++i) \
    for (size_t j = 0; j < (c); ++j)

#define CHECK_SIZE(mat)                                          \
  do {                                                           \
    if (!((mat).nRows() == nRows() && (mat).nCols() == nCols())) \
      throw incompatible_size();                                 \
  } while (0)

#define CHECK_ROW(r)              \
  do {                            \
    if ((r) < 0 || (r) > nRows()) \
      throw row_outbound();       \
  } while (0)

#define CHECK_COLUMN(c)           \
  do {                            \
    if ((c) < 0 || (c) > nCols()) \
      throw col_outbound();       \
  } while (0)
#define MATCH_ROW(r)         \
  do {                       \
    if ((r) != nRows())      \
      throw bad_row_match(); \
  } while (0)
#define MATCH_COLUMN(c)      \
  do {                       \
    if ((c) != nCols())      \
      throw bad_col_match(); \
  } while (0)
#define VALIDATE_INDEX(r, c) \
  do {                       \
    CHECK_ROW(r);            \
    CHECK_COLUMN(c);         \
  } while (0)

template <typename T>
class Matrix {
 private:
  size_t mRows;
  size_t mCols;
  size_t mSize;
  std::vector<T> mData;

 public:
  size_t nCols() const { return mCols; }
  size_t nRows() const { return mRows; }
  size_t size() const { return mSize; }
  const std::vector<T>& data() const { return mData; }
  void updateSize() { mSize = mRows * mCols; }

  T& operator()(size_t row, size_t col) {
    try {
      validateIndexes(row, col);

      return mData[row * nCols() + col];
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;

      throw;
    }
  }

  const T& operator()(size_t row, size_t col) const {
    try {
      validateIndexes(row, col);
      return mData[row * nCols() + col];
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;

      throw;
    }
  }

  // Constructors and destructor...
  Matrix() : mRows(0), mCols(0), mSize(0), mData(0) {}

  Matrix(size_t dimension)
      : mRows(dimension),
        mCols(dimension),
        mSize(dimension * dimension),
        mData(dimension * dimension, 0) {}

  Matrix(size_t rows, size_t cols)
      : mRows(rows), mCols(cols), mSize(rows * cols), mData(rows * cols, 0) {}
  //copy constructor
  Matrix(const Matrix& matrix)
      : mRows(matrix.mRows),
        mCols(matrix.mCols),
        mSize(matrix.mSize),
        mData(matrix.mData) {}

  //constructor for converting one type to other
  template <typename U>
  Matrix(const Matrix<U>& other)
      : mRows(other.nRows()),
        mCols(other.nCols()),
        mSize(mRows * mCols),
        mData(mSize) {
    try {
      for (size_t i = 0; i < mRows; ++i) {
        for (size_t j = 0; j < mCols; ++j) {
          (*this)(i, j) =
              static_cast<T>(other(i, j));  // Convert each element to type T
        }
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }

  //constructor from a initializer list
  Matrix(std::initializer_list<std::initializer_list<T>> list)
      : mRows(list.size()), mCols(list.begin()->size()), mSize(mRows * mCols) {
    try {
      for (const auto& row : list) {
        if (row.size() != mCols) {
          throw bad_rows();
        }
        mData.insert(mData.end(), row.begin(), row.end());
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }
  //constructor from a initializer list  and given adequate row and column
  Matrix(size_t rows, size_t cols, std::initializer_list<T> initList)
      : mRows(rows), mCols(cols), mSize(rows * cols), mData(initList) {
    try {
      if (initList.size() != mSize) {
        throw bad_size();
      }

    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }
  //constructor to create matrix from std::vector
  Matrix(size_t rows, size_t cols, const std::vector<T>& data)
      : mRows(rows), mCols(cols), mSize(rows * cols), mData(data) {
    try {
      if (data.size() != mSize) {
        throw bad_size();
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }
  //constructor to create matrix from arrays
  template <std::size_t N>
  Matrix(size_t rows, size_t cols, T (&data)[N])
      : mRows(rows), mCols(cols), mSize(rows * cols), mData(data, data + N) {
    try {
      if (N != mSize) {
        throw bad_size();
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }
  // Constructor for row or column vector
  Matrix(const std::vector<T>& data, VectorType vectorType) : mData(data) {
    try {
      if (vectorType == VectorType::RowVector) {
        mRows = 1;
        mCols = data.size();
      } else if (vectorType == VectorType::ColumnVector) {
        mRows = data.size();
        mCols = 1;
      } else {
        throw bad_argument();
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }
  //constructor to create matrix from list of vectors
  Matrix(const std::vector<Matrix<T>>& vectors, Orientation orientation) {
    try {
      if (orientation == Orientation::Row) {
        mRows = vectors.size();
        mCols = mRows > 0 ? vectors[0].nCols() : 0;
        // Check if all vectors are row vectors
        for (const auto& vec : vectors) {
          if (vec.nCols() != mCols)
            throw bad_size();
          if (vec.nRows() != 1)
            throw too_many_rows();
        }
      } else if (orientation == Orientation::Column) {
        mCols = vectors.size();
        mRows = mCols > 0 ? vectors[0].nRows() : 0;
        // Check if all vectors are column vectors
        for (const auto& vec : vectors) {
          if (vec.nRows() != mRows)
            throw bad_size();
          if (vec.nCols() != 1)
            throw too_many_cols();
        }
      } else {
        throw bad_argument();
      }

      // Combine the vectors into the matrix
      for (const auto& vec : vectors) {
        for (const auto& item : vec.data()) {
          mData.push_back(item);
        }
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;

      throw;  // Re-throw the caught exception by reference
    }
  }

  // Destructor
  ~Matrix() {}
  // Scalar operators
  friend Matrix operator+(const Matrix& m, double scalar) {
    Matrix result(m.nRows(), m.nCols());
#pragma omp parallel for collapse(2)
    for_ij(m.nRows(), m.nCols()) result(i, j) = m(i, j) + scalar;
    return result;
  }

  friend Matrix operator+(double scalar, const Matrix& m) { return m + scalar; }

  friend Matrix operator-(const Matrix& m, double scalar) {
    return m + (-scalar);
  }

  friend Matrix operator-(double scalar, const Matrix& m) {
    return -(m - scalar);
  }

  friend Matrix operator*(const Matrix& m, double scalar) {
    Matrix result(m.nRows(), m.nCols());
#pragma omp parallel for collapse(2)
    for_ij(m.nRows(), m.nCols()) result(i, j) = m(i, j) * scalar;
    return result;
  }

  friend Matrix operator*(double scalar, const Matrix& m) { return m * scalar; }

  friend Matrix operator/(const Matrix& m, double scalar) {
    // Check for division by zero
    if (scalar == 0.0) {
      throw std::invalid_argument("Division by zero");
    }

    Matrix result(m.nRows(), m.nCols());

// Parallelize the loop using OpenMP with collapse(2)
#pragma omp parallel for collapse(2)
    for_ij(m.nRows(), m.nCols()) {
      result(i, j) = m(i, j) / scalar;
    }

    return result;
  }

  friend Matrix operator/(double scalar, const Matrix& m) {
    // division is not commutative, so a new method is implemented
    Matrix result(m.nRows(), m.nCols());
#pragma omp parallel for collapse(2)
    for_ij(m.nRows(), m.nCols()) result(i, j) = scalar / m(i, j);
    return result;
  }
  Matrix& operator+=(T scalar);
  //   Matrix& operator+=(T scalar) {
  //     // Perform addition element-wise with the scalar value
  //     // Parallelization pragma
  //     #pragma omp parallel for collapse(2)
  //     for_ij(nRows(), nCols()) {
  //         (*this)(i, j) += scalar;
  //     }
  //     return *this;
  // }

  //   Matrix operator+=(double scalar) {
  // #pragma omp parallel for
  //     for_index(i, mData.size()) mData[i] += scalar;
  //     return *this;
  //   }

  Matrix operator-=(double scalar) {
#pragma omp parallel for
    for_index(i, mData.size()) mData[i] -= scalar;

    return *this;
  }

  Matrix operator*=(double scalar) {
#pragma omp parallel for
    for_index(i, mData.size()) mData[i] *= scalar;
    return *this;
  }

  Matrix& operator/=(double scalar) {
    // Check for division by zero
    if (scalar == 0.0) {
      throw std::invalid_argument("Division by zero");
    }

// Parallelize the loop using OpenMP
#pragma omp parallel for
    for_index(i, size()) {
      mData[i] /= scalar;
    }

    return *this;
  }

  Matrix<T> operator==(const T& scalar) {
    Matrix<T> result(nRows(), nCols());
#pragma omp parallel for collapse(2)
    for_ij(nRows(), nCols()) result(i, j) = operator()(i, j) == scalar;
    return result;
  }

  Matrix operator!=(const double& scalar) { return -((*this == scalar) - 1); }
  // endregion scalar operator

  //-------------------------- region Matrix operators--------------------------------------//

  template <typename T2>
  auto operator+(const Matrix<T2>& m) const {
    // Determine the superior type using std::common_type
    using result_type = typename std::common_type<T, T2>::type;

    // Perform size checking using the macro
    CHECK_SIZE(m);

    Matrix<result_type> result(nRows(), nCols());
// Parallelization pragma
#pragma omp parallel for collapse(2)
    for_ij(nRows(), nCols()) {
      result(i, j) = static_cast<result_type>((*this)(i, j)) +
                     static_cast<result_type>(m(i, j));
    }
    return result;
  }

  template <typename T2>
  auto operator-(const Matrix<T2>& m) const {
    // Determine the superior type using std::common_type
    using result_type = typename std::common_type<T, T2>::type;

    // Perform size checking using the macro
    CHECK_SIZE(m);

    Matrix<result_type> result(nRows(), nCols());
// Parallelization pragma
#pragma omp parallel for collapse(2)
    for_ij(nRows(), nCols()) {
      result(i, j) = static_cast<result_type>((*this)(i, j)) -
                     static_cast<result_type>(m(i, j));
    }
    return result;
  }

  template <typename T2>
  auto operator*(const Matrix<T2>& m2) const {
    // Determine the superior type using std::common_type
    using result_type = typename std::common_type<T, T2>::type;

    Matrix<result_type> result(nRows(), m2.nCols());
// Parallelization pragma
#pragma omp parallel for if (result.nRows() * result.nCols() > 250)
    for (size_t i = 0; i < nRows(); ++i) {
      for (size_t j = 0; j < m2.nCols(); ++j) {
        for (size_t k = 0; k < nCols(); ++k) {
          result(i, j) += static_cast<result_type>((*this)(i, k)) *
                          static_cast<result_type>(m2(k, j));
        }
      }
    }
    return result;
  }

  template <typename T2>
  Matrix& operator+=(const Matrix<T2>& m) {
    // Ensure matrices have the same size
    if (nRows() != m.nRows() || nCols() != m.nCols())
      throw std::invalid_argument("Matrix sizes do not match");

    // Determine the result type of the addition using std::common_type
    using result_type = typename std::common_type<T, T2>::type;
    // Promote to Matrix<float> if T is int
    // if constexpr (std::is_same_v<T, int>) {
    //  // std::cout<<"hello"<<std::endl;
    //     *this = static_cast<Matrix<float>>(*this);

    // }
    Matrix<float> temp(*this);
    // Perform addition element-wise
    for (size_t i = 0; i < nRows(); ++i) {
      for (size_t j = 0; j < nCols(); ++j) {

        // Convert both the left-hand side and right-hand side to the result type
        temp(i, j) += static_cast<float>(m(i, j));
      }
    }
    *this = temp;
    return (*this);
  }

  // Define the template function for the addition assignment operator
  //  template<typename T2>
  //   Matrix& operator+=(const Matrix<T2>& m) {
  //       // Perform size checking using the macro
  //       CHECK_SIZE(m);

  //       // Determine the superior type using std::common_type
  //       using result_type = typename std::common_type<T, T2>::type;

  //       // Perform addition element-wise
  //       // Parallelization pragma
  //       #pragma omp parallel for collapse(2)
  //       for_ij(nRows(), nCols()) {
  //           (*this)(i, j) += static_cast<T>(static_cast<result_type>((*this)(i, j)) + static_cast<result_type>(m(i, j)));
  //       }

  //       return *this;
  //   }

  // template<typename T2>
  //   Matrix& operator+=(const Matrix<T2>& m);

  //   Matrix& operator+=(const Matrix& other) {
  //     CHECK_SIZE(other);
  // #pragma omp parallel for collapse(2)
  //     for_ij(other.nRows(), other.nCols()) operator()(i, j) += other(i, j);
  //     return *this;
  //   }

  Matrix& operator-=(const Matrix& other) {
    CHECK_SIZE(other);

#pragma omp parallel for collapse(2)
    for_ij(other.nRows(), other.nCols()) operator()(i, j) -= other(i, j);
    return *this;
  }

  Matrix& operator*=(const Matrix& other) {
    if (nCols() != other.nRows())
      throw incompatible_size();

    Matrix result(nRows(), other.nCols());

#pragma omp parallel for collapse(2)
    // two loops iterate through every cell of the new matrix
    for (size_t i = 0; i < result.nRows(); i++) {
      for (size_t j = 0; j < result.nCols(); j++) {
        // here we calculate the value of a single cell in our new matrix
        result(i, j) = 0;
        for (size_t k = 0; k < nCols(); k++)
          result(i, j) += operator()(i, k) * other(k, j);
      }
    }

    mRows = result.nRows();
    mCols = result.nCols();
    mData = result.mData;
    return *this;
  }

  bool operator==(const Matrix& other) {
    if (mData.size() != other.mData.size() || mRows != other.mRows ||
        mCols != other.mCols)
      return false;

    for (int k = 0; k < mData.size(); k++) {
      if (mData[k] != other.mData[k])
        return false;
    }

    return true;
  }

  bool operator!=(const Matrix& other) { return !(*this == other); }

  Matrix slice(size_t start, size_t end, Orientation sliceType) const {
    try {
      if (sliceType == Orientation ::Row) {
        if (start >= mRows || end > mRows || start > end) {
          throw std::out_of_range("Invalid row slice indices");
        }
        size_t sliceRows = end - start;
        Matrix result(sliceRows, mCols);
        for (size_t i = start; i < end; ++i) {
          for (size_t j = 0; j < mCols; ++j) {
            result(i - start, j) = (*this)(i, j);
          }
        }
        return result;
      } else if (sliceType == Orientation ::Column) {
        if (start >= mCols || end > mCols || start > end) {
          throw std::out_of_range("Invalid column slice indices");
        }
        size_t sliceCols = end - start;
        Matrix result(mRows, sliceCols);
        for (size_t i = 0; i < mRows; ++i) {
          for (size_t j = start; j < end; ++j) {
            result(i, j - start) = (*this)(i, j);
          }
        }
        return result;
      } else {
        throw std::invalid_argument("Invalid slice type");
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }

  Matrix operator[](
      std::pair<std::pair<size_t, size_t>, Orientation> indices) const {
    size_t start = indices.first.first;
    size_t end = indices.first.second;
    Orientation sliceType = indices.second;

    if (sliceType == Orientation::Row) {
      return slice(start, end, Orientation::Row);
    } else if (sliceType == Orientation::Column) {
      return slice(start, end, Orientation::Column);
    } else {
      throw std::invalid_argument("Invalid slice type");
    }
  }

  Matrix operator[](std::pair<size_t, Orientation> index) const {
    size_t idx = index.first;
    Orientation sliceType = index.second;

    if (sliceType == Orientation::Row) {
      return slice(idx, idx + 1, Orientation::Row);
    } else if (sliceType == Orientation::Column) {
      return slice(idx, idx + 1, Orientation::Column);
    } else {
      throw std::invalid_argument("Invalid slice type");
    }
  }

  // endregion matrix operator

  class iterator {
   private:
    size_t mRow;
    size_t mCol;
    Matrix<T>* mMatrixPtr;

   public:
    size_t getrow() { return mRow; }

    size_t getcol() { return mCol; }

    iterator(size_t row, size_t col, Matrix<T>* matrixPtr)
        : mRow(row), mCol(col), mMatrixPtr(matrixPtr) {}

    T& operator*() const {
      if (mRow >= mMatrixPtr->nRows() || mCol >= mMatrixPtr->nCols()) {
        throw out_of_range();
      }
      return (*mMatrixPtr)(mRow, mCol);
    }

    bool operator==(const iterator& other) const {
      return (mRow == other.mRow && mCol == other.mCol);
    }

    bool operator!=(const iterator& other) const { return !(*this == other); }

    // Comparison operator <
    bool operator<(const iterator& other) const {
      // Compare row and column indices lexicographically
      if (mRow != other.mRow)
        return mRow < other.mRow;
      else
        return mCol < other.mCol;
    }

    // Comparison operator <=
    bool operator<=(const iterator& other) const {
      // Compare row and column indices lexicographically
      if (mRow != other.mRow)
        return mRow < other.mRow;
      else
        return mCol <= other.mCol;
    }

    iterator& operator++() {
      ++mCol;
      if (mCol == mMatrixPtr->nCols()) {
        ++mRow;
        mCol = 0;
      }
      return *this;
    }

    iterator operator++(int) {
      iterator temp = *this;
      ++(*this);
      return temp;
    }

    iterator& operator--() {
      if (mCol == 0) {
        --mRow;
        mCol = mMatrixPtr->nCols() - 1;
      } else {
        --mCol;
      }
      return *this;
    }

    iterator operator--(int) {
      iterator temp = *this;
      --(*this);
      return temp;
    }

    // Advance function
    void advance(int n) {
      if (n > 0) {
        for (int i = 0; i < n; ++i) {
          ++(*this);
        }
      } else if (n < 0) {
        for (int i = 0; i < -n; ++i) {
          --(*this);
        }
      }
    }

    // Distance to function
    int distance_to(const iterator& other) const {
      if (mMatrixPtr != other.mMatrixPtr) {
        throw std::invalid_argument("Iterators belong to different matrices.");
      }

      int currentIndex = static_cast<int>(mRow * mMatrixPtr->nCols() + mCol);
      int otherIndex =
          static_cast<int>(other.mRow * mMatrixPtr->nCols() + other.mCol);
      return otherIndex - currentIndex;
    }
  };

  iterator begin() { return iterator(0, 0, this); }

  iterator end() { return iterator(mRows, 0, this); }

  typename std::vector<T>::iterator row_begin(size_t row) {
    if (row >= mRows) {
      throw row_outbound();
    }
    return mData.begin() + (row * mCols);
  }

  // Get iterator to the end of a row
  typename std::vector<T>::iterator row_end(size_t row) {
    if (row >= mRows) {
      throw row_outbound();
    }
    return mData.begin() + ((row + 1) * mCols);
  }

  // Get iterator to the beginning of a column
  typename std::vector<T>::iterator col_begin(size_t col) {
    if (col >= mCols) {
      throw col_outbound();
    }
    return mData.begin() + col;
  }

  // Get iterator to the end of a column
  typename std::vector<T>::iterator col_end(size_t col) {
    if (col >= mCols) {
      throw col_outbound();
    }
    return mData.end() - (mCols - col);
  }

  void validateIndexes(size_t row, size_t col) const {

    if (row >= mRows || col >= mCols)
      throw std::out_of_range("Index out of bounds");
  }

  void validateInsertIndexes(size_t index, size_t insert_dim, size_t other_dim,
                             Orientation orientation) const {
    if (orientation == Orientation::Row) {
      if (insert_dim != 1) {
        throw too_many_rows();
      }
      if (other_dim != mCols) {
        throw bad_rows();
      }
      if (index > mRows) {
        throw row_outbound();
      }
    } else if (orientation == Orientation::Column) {
      if (insert_dim != 1) {
        throw too_many_cols();
      }
      if (other_dim != mRows) {
        throw bad_cols();
      }
      if (index > mCols) {
        throw col_outbound();
      }
    }
  }

  // Static member functions for creating special matrices
  static Matrix<T> identity(size_t size) {
    Matrix<T> identityMatrix(size, size);
    for (size_t i = 0; i < size; ++i) {
      for (size_t j = 0; j < size; ++j) {
        identityMatrix(i, j) = (i == j) ? static_cast<T>(1) : static_cast<T>(0);
      }
    }
    return identityMatrix;
  }

  static Matrix<T> zeros(size_t rows, size_t cols) {
    return Matrix<T>(rows, cols);
  }

  static Matrix<T> ones(size_t rows, size_t cols) {
    Matrix<T> result(rows, cols);
    for (auto& element : result.mData) {
      element = static_cast<T>(1);
    }
    return result;
  }

  static Matrix<T> diagonal(const std::vector<T>& elements) {
    size_t dimension = elements.size();
    Matrix<T> result(dimension, dimension);
    for (size_t i = 0; i < dimension; ++i) {
      result(i, i) = elements[i];
    }
    return result;
  }
  // Getter functions...

  Matrix operator-() {
    Matrix result(this->mRows, this->mCols);

#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < mCols; i++) {
      for (size_t j = 0; j < mRows; j++) {
        result(i, j) = -operator()(i, j);
      }
    }

    return result;
  }
  // Getter function
  T get_element(size_t row, size_t col) const {
    try {
      validateIndexes(row, col);
      return mData[row * nCols() + col];
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;

      throw;
    }
  }
  // Setter function
  void set_element(size_t row, size_t col, const T& value) {
    try {
      validateIndexes(row, col);
      mData[row * nCols() + col] = value;
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }
  template <typename U>
  Matrix<U> augmentedMatrix(const Matrix<U>& other) const {
    MATCH_ROW(other.nRows());
    Matrix<U> augmented(this->nRows(), this->nCols() + other.nCols());

    // Copy the original matrix into the left half of the augmented matrix
    for (size_t i = 0; i < this->nRows(); ++i) {
      for (size_t j = 0; j < this->nCols(); ++j) {
        augmented(i, j) = this->mData[i * this->nCols() + j];
      }
    }

    // Copy the other matrix into the right half of the augmented matrix
    for (size_t i = 0; i < this->nRows(); ++i) {
      for (size_t j = 0; j < other.nCols(); ++j) {
        augmented(i, this->nCols() + j) = other(i, j);
      }
    }

    return augmented;
  }

  template <typename U>
  Matrix<U> augmentedMatrix(const std::vector<U>& vector) const {
    MATCH_ROW(vector.size());

    Matrix<U> augmented(this->nRows(), this->nCols() + 1);

    // Copy the original matrix into the left half of the augmented matrix
    for (size_t i = 0; i < this->nRows(); ++i) {
      for (size_t j = 0; j < this->nCols(); ++j) {
        augmented(i, j) = this->mData[i * this->nCols() + j];
      }
    }

    // Copy the vector into the rightmost column of the augmented matrix
    for (size_t i = 0; i < this->nRows(); ++i) {
      augmented(i, this->nCols()) = vector[i];
    }

    return augmented;
  }

  // region basic matrix  boolean functions
  bool isEmpty() { return mCols == 0 and mRows == 0; }

  bool isSquare() const { return nRows() == nCols(); }
  bool isColumn() const { return mCols == 1; }

  bool isRow() const { return mRows == 1; }

  bool contains(T value) {
    return std::find(mData.begin(), mData.end(), value) != mData.end();
  }
  //   Matrix diagonal() {
  //     if (!isSquare()) {
  //       throw std::runtime_error("Can't get the diagonal, not a square matrix");
  //     }
  //     Matrix result(mRows, 1);

  // #pragma omp parallel
  //     for (size_t i = 0; i < mRows; i++)
  //       result(i, 0) = operator()(i, i);

  //     return result;
  //   }

  std::tuple<size_t, size_t, T> min() const {
    T minVal = std::numeric_limits<T>::max();
    // T minVal = mData[0];
    size_t minPos = 0;

    for (size_t i = 0; i < mRows; ++i) {
      for (size_t j = 0; j < mCols; ++j) {
        if (minVal > mData[i * mCols + j]) {
          minPos = i * mCols + j;
          minVal = mData[minPos];
        }
      }
    }

    size_t r = minPos / mCols;  // Calculate row index
    size_t c = minPos % mCols;  // Calculate column index

    return std::make_tuple(r, c, minVal);
  }

  T min_element() const {
    return *std::min_element(std::begin(mData), std::end(mData));
  }

  T max() const {
    return *std::max_element(std::begin(mData), std::end(mData));
  }

  T sum() const {
    T sum_of_elems = 0;
    for (T n : mData)
      sum_of_elems += n;

    return sum_of_elems;
  }
  // end region functors

  // region functions (row and column operation)
  Matrix getRow(size_t index) {
    if (index >= mRows)
      throw std::invalid_argument("Row index out of bounds");

    Matrix result(1, mCols);
#pragma omp parallel for
    for (size_t i = 0; i < mCols; i++)
      result(0, i) = operator()(index, i);

    return result;
  }

  Matrix getColumn(size_t index) {
    if (index >= mCols)
      throw std::invalid_argument("Column index out of bounds");

    Matrix result(mRows, 1);
#pragma omp parallel for
    for (size_t i = 0; i < mRows; i++)
      result(i, 0) = operator()(i, index);

    return result;
  }

  template <typename InputIterator, typename OutputIterator>
  void insertSubmatrix(InputIterator srcBegin, InputIterator srcEnd,
                       OutputIterator destPos, size_t numRows, size_t numCols,
                       size_t position, size_t offset) {
    // Copy the values from source to destination with the appropriate offset
    for (size_t i = 0; i < numRows; ++i) {
      std::copy(srcBegin + i * numCols, srcBegin + (i + 1) * numCols,
                destPos + (position + i) * offset);
    }
  }

  void addRow(const Matrix& values, size_t position) {
    try {
      if (isEmpty()) {
        mRows = 1;  // Initialize to 1 row initially
        mCols = values.nCols();
        mData = values.mData;
        return;
      }

      validateInsertIndexes(position, values.nRows(), values.nCols(),
                            Orientation::Row);
      // Prepare a new data vector with space for the additional new row
      auto rowBegin = mData.begin() + position * mCols;
      // Insert the new row at the specified position
      mData.insert(rowBegin, values.mData.begin(), values.mData.end());

      mRows += 1;
      updateSize();
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return;
    }
  }

  void addRow(std::initializer_list<T> rowValues, size_t position) {
    // Validate input
    try {
      if (isEmpty()) {
        mRows = 1;  // Initialize to 1 row initially
        mCols = rowValues.size();
        mData.insert(mData.end(), rowValues.begin(), rowValues.end());
        return;
      } else {
        validateInsertIndexes(position, 1, mCols, Orientation::Row);
        auto rowBegin = mData.begin() + position * mCols;
        mData.insert(rowBegin, rowValues.begin(), rowValues.end());
        mRows += 1;
        updateSize();
      }
    } catch (const std::invalid_argument& e) {
      std::cerr << "Error: " << e.what() << std::endl;

      return;
    }
  }

  void addColumn(const Matrix& values, size_t position) {
    try {
      if (isEmpty()) {
        mRows = values.nRows();
        mCols = 1;
        mData = values.mData;
        return;
      }

      validateInsertIndexes(position, values.nCols(), values.nRows(),
                            Orientation::Column);

      // Insert the new column at the specified position
      for (size_t i = 0; i < mRows; ++i) {
        mData.insert(mData.begin() + (i * (mCols + 1)) + position,
                     values(i, 0));
      }

      mCols += 1;
      updateSize();
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return;
    }
  }

  void addColumn(std::initializer_list<T> colValues, size_t position) {
    // Validate input
    try {
      if (isEmpty()) {
        mRows = colValues.size();
        mCols = 1;  // Initialize to 1 column initially
        for (const auto& value : colValues) {
          mData.push_back(value);
        }
        return;
      } else {
        validateInsertIndexes(position, 1, mRows, Orientation::Column);
        for (size_t i = 0; i < mRows; ++i) {
          mData.insert(mData.begin() + (i * (mCols + 1)) + position,
                       *(colValues.begin() + i));
        }
        mCols += 1;
        updateSize();
      }
    } catch (const std::invalid_argument& e) {
      std::cerr << "Error: " << e.what() << std::endl;

      return;
    }
  }

  void addRow(Matrix values) { addRow(values, mRows); }

  void addRow(std::initializer_list<T> rowValues) { addRow(rowValues, mRows); }

  void addColumn(Matrix values) { addColumn(values, mCols); }

  void addColumn(std::initializer_list<T> colValues) {
    addColumn(colValues, mCols);
  }

  Matrix<T>& setRow(Matrix<T> rowVector, size_t rowPos) {
    MATCH_COLUMN(rowVector.nCols());
    CHECK_ROW(rowPos);

    if (rowVector.nRows() != 1)
      throw too_many_rows();

    for (size_t col = 0; col < nCols(); col++)
      operator()(rowPos, col) = rowVector(0, col);

    return *this;
  }

  Matrix<T>& setColumn(Matrix<T> columnVector, size_t colPos) {
    MATCH_ROW(columnVector.nRows());
    CHECK_COLUMN(colPos);
    if (columnVector.nCols() > 1)
      throw too_many_cols();

    for (size_t row = 0; row < nRows(); row++)
      operator()(row, colPos) = columnVector(row, 0);
    return *this;
  }

  Matrix<T>& reshape(size_t rows, size_t cols) {
    CHECK_SIZE(rows * cols);
    mRows = rows;
    mCols = cols;
    updateSize();
    return *this;
  }

  Matrix<T>& swapRows(size_t r1, size_t r2) {
    if (std::max(r1, r2) >= nRows())
      throw row_outbound();

    for (size_t col = 0; col < nCols(); col++) {
      std::swap(operator()(r1, col), operator()(r2, col));
    }

    return *this;
  }

  Matrix<T>& swapColumns(size_t c1, size_t c2) {
    if (std::max(c1, c2) >= nCols())
      throw col_outbound();

    for (size_t row = 0; row < nRows(); row++) {
      std::swap(operator()(row, c1), operator()(row, c2));
    }

    return *this;
  }

  Matrix<T>& removeRow(size_t position) {
    CHECK_ROW(position);
    mData.erase(mData.begin() + position * nCols(),
                mData.begin() + (position + 1) * nCols());
    --mRows;
    updateSize();
    return *this;
  }

  Matrix<T>& removeColumn(size_t position) {
    CHECK_COLUMN(position);
    for (auto it = mData.begin() + position; it != mData.end();) {
      it = mData.erase(it);
      it += nCols() - 1;  // Move the iterator to the next column
    }
    --mCols;
    updateSize();
    return *this;
  }

  Matrix filter(const Matrix<T>& bin, bool columns = false) {
    size_t dimension = columns ? mCols : mRows;
    Matrix result;

    try {
      if (bin.nCols() != 1)
        throw too_many_cols();
      if (bin.nRows() > dimension)
        throw bad_size();

      for (size_t i = 0; i < bin.nRows(); i++) {
        if (bin(i, 0)) {
          if (columns)
            result.addColumn(getColumn(i));
          else
            result.addRow(getRow(i));
        }
      }
    } catch (const std::invalid_argument& e) {
      std::cerr << "Error: " << e.what() << std::endl;
    }

    return result;
  }

  Matrix getRows(const Matrix<T> bin) { return filter(bin); }

  Matrix getColumns(const Matrix<T> bin) { return filter(bin, true); }

  // region elementary matrix

  template <typename U>
  Matrix<U> fill(U value) const {
    Matrix<U> filledMatrix(nRows(), nCols());
    for_ij(nRows(), nCols()) filledMatrix(i, j) = value;
    return filledMatrix;
  }

  Matrix<T> to_ones() { return fill(T(1)); }

  Matrix<T> to_zeros() { return fill(T(0)); }

  Matrix<T> to_diagonal(T value) {
    size_t size = std::min(nRows(), nCols());
    Matrix<T> result(size, size);
    for_index(i, size) result(i, i) = value;
    return result;
  }

  Matrix<T> to_identity() { return to_diagonal(1); }

  Matrix E_SwapOp(size_t r0, size_t r1) {
    try {
      if (std::max(r0, r1) >= mRows) {
        throw std::out_of_range("row index exceeds matrix size");
      }

      Matrix<T> swapOp = to_identity();

      // Swap rows r0 and r1
      for (size_t j = 0; j < mCols; ++j) {
        std::swap(swapOp(r0, j), swapOp(r1, j));
      }

      return swapOp;
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }

  Matrix<T> E_MultiplyOp(T factor, size_t r) {
    try {
      if (r >= nRows()) {
        throw std::out_of_range("row index exceeds matrix size");
      }

      Matrix<T> mulOp = to_identity();

      mulOp(r, r) = factor;

      return mulOp;
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }

  Matrix<T> E_subtract_product_row_op(T factor, size_t r0, size_t r1) {
    try {
      if (std::max(r0, r1) >= nRows()) {
        throw std::out_of_range("Row index exceeds matrix size");
      }

      Matrix<T> productOp = to_identity();

      // Subtract factor times row r0 from row r1
      for (size_t j = 0; j < nCols(); ++j) {
        productOp(r1, j) -= factor * productOp(r0, j);
      }

      return productOp;
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }

  template <typename T>
  void normalize_pivot_row(Matrix<T>& matrix, const T& pivotElement,
                           size_t processingRow,
                           size_t processingColumn) const {
    if (pivotElement == 0) {
      return;  // Skip normalization if pivot element is zero
    }

    // Make pivot row 1 and zero out other entries in the pivot column
    for (size_t column = processingColumn; column < matrix.nCols(); ++column) {
      matrix(processingRow, column) /= pivotElement;
      // Check if the value is close to zero and set it to zero if so
      if (std::abs(matrix(processingRow, column)) < 1e-10) {
        matrix(processingRow, column) = 0;
      }
    }
  }

  template <typename T>
  void eliminate_rows(Matrix<T>& matrix, const T& pivotElement,
                      size_t processingRow, size_t processingColumn,
                      bool reduced_form = false) const {
    if (reduced_form) {
      // Eliminate all rows except processing row above or below pivot for reduced form
      for (size_t row = 0; row < matrix.nRows(); ++row) {
        T localPivotElement = matrix(row, processingColumn);

        T scale_factor = localPivotElement;
        if (row != processingRow && localPivotElement != 0) {
          for (size_t column = processingColumn; column < matrix.nCols();
               ++column) {
            matrix(row, column) -= matrix(processingRow, column) * scale_factor;
            // Check if the value is close to zero and set it to zero if so
            if (EQUAL(std::abs(matrix(row, column)), 0)) {
              matrix(row, column) = 0;
            }
          }
        }
      }
    } else {
      // Eliminate rows below the processing row
      if (pivotElement == 0) {
        return;  // Skip elimination if pivot element is zero
      }
      for (size_t row = processingRow + 1; row < matrix.nRows(); ++row) {
        T scale_factor = matrix(row, processingColumn) / pivotElement;
        for (size_t column = processingColumn; column < matrix.nCols();
             ++column) {
          matrix(row, column) -= matrix(processingRow, column) * scale_factor;
          // Check if the value is close to zero and set it to zero if so
          if (EQUAL(std::abs(matrix(row, column)), 0)) {
            matrix(row, column) = 0;
          }
        }
      }
    }
  }

  std::pair<Matrix<T>, int> echelonizeWithSign(
      bool reduced_row_echelon = false) const {
    Matrix<T> echelonMatrix(*this);  // Create a copy of the current matrix
    int sign = 1;                    // Initialize sign
    size_t processingRow = 0;
    for (size_t column = 0; column < nCols(); ++column) {
      // Find the pivot element in the current column
      size_t pivotRow = processingRow;
      while (pivotRow < nRows() && echelonMatrix(pivotRow, column) == 0) {
        ++pivotRow;
      }

      if (pivotRow == nRows()) {
        continue;  // No pivot found, move to the next column
      }

      if (pivotRow != processingRow) {
        // Swap rows to bring the pivot element to the current processing row
        echelonMatrix.swapRows(processingRow, pivotRow);
        sign *= -1;  // Update sign due to row swap
      }
      T pivotElement = echelonMatrix(pivotRow, column);
      eliminate_rows(echelonMatrix, pivotElement, processingRow, column);

      // Make pivot row 1 and zero out other entries in the pivot column
      if (reduced_row_echelon && pivotElement != 0) {
        normalize_pivot_row(echelonMatrix, pivotElement, processingRow, column);
        eliminate_rows(echelonMatrix, pivotElement, processingRow, column,
                       true);
      }

      processingRow++;
    }

    return std::make_pair(echelonMatrix, sign);
  }

  Matrix<T> rref2() {
    Matrix<T> result(*this);

    size_t ProcessColumn = 0;
    const size_t rowCount = result.nRows();
    const size_t colCount = result.nCols();

    for (size_t row = 0; row < rowCount; row++) {
      if (ProcessColumn >= colCount) {
        return result;
      }

      // Find the pivot row
      size_t pivotRow = row;
      while (EQUAL(result(pivotRow, ProcessColumn), 0)) {
        pivotRow++;
        if (pivotRow == rowCount) {
          pivotRow = row;
          ProcessColumn++;
          if (ProcessColumn == colCount) {
            return result;
          }
        }
      }
      result.swapRows(pivotRow, row);
      T pivotValue = result(row, ProcessColumn);
      // Scale the pivot row to have a leading 1
      if (!EQUAL(pivotValue, 0)) {

        auto tmp = result.getRow(row) /= pivotValue;

        result.setRow(tmp, row);
      }

      // Eliminate all other entries in the current column
      for (size_t other_row = 0; other_row < rowCount; other_row++) {
        if (other_row != row) {
          T factor = result(other_row, ProcessColumn);
          for (size_t column = 0; column < colCount; column++) {
            result(other_row, column) -= factor * result(row, column);
          }
        }
      }

      ProcessColumn++;
    }

    return result;
  }

  Matrix<T> rref3() {
    Matrix<T> result(*this);

    size_t ProcessColumn = 0;
    const size_t rowCount = result.nRows();
    const size_t colCount = result.nCols();

    for (size_t row = 0; row < rowCount; row++) {
      if (ProcessColumn >= colCount) {
        return result;
      }

      // Find the pivot row
      size_t pivotRow = row;
      while (pivotRow < rowCount && EQUAL(result(pivotRow, ProcessColumn), 0)) {
        pivotRow++;
      }
      if (pivotRow == rowCount) {
        ProcessColumn++;
        continue;
      }

      // Scale the pivot row to have a leading 1
      T pivotValue = result(pivotRow, ProcessColumn);
      if (!EQUAL(pivotValue, 0)) {
        auto tmp = result.getRow(pivotRow) /= pivotValue;
        result.setRow(tmp, pivotRow);
      }

      // Eliminate all other entries in the current column
      for (size_t other_row = 0; other_row < rowCount; other_row++) {
        if (other_row != pivotRow) {
          T factor = result(other_row, ProcessColumn);
          for (size_t column = 0; column < colCount; column++) {
            result(other_row, column) -= factor * result(pivotRow, column);
          }
        }
      }

      ProcessColumn++;
    }

    return result;
  }

  Matrix<T> rowEchelon() const {

    auto [echelonMatrix, sign] = echelonizeWithSign();
    return echelonMatrix;  // Call echelon with default argument
  }

  Matrix<T> rref() {
    auto [echelonMatrix, sign] = echelonizeWithSign(true);
    return echelonMatrix;  // Call echelon with make_pivot_one set to true
  }

  struct Decomposition_LU {
    Matrix<T> L;
    Matrix<T> U;
  };

  struct Decomposition_LDU {
    Matrix<T> L;
    Matrix<T> D;
    Matrix<T> U;
  };

  Decomposition_LU to_LU() const {
    Matrix<T> A(*this);  // Make a copy of the matrix

    size_t n = A.nRows();
    Matrix<T> L(n, n);             // Lower triangular matrix
    Matrix<T> U = A.rowEchelon();  // Upper triangular matrix

    // Compute lower triangular matrix L
    for (size_t i = 0; i < n; ++i) {
      L(i, i) = 1;  // Diagonal elements of L are set to 1
      for (size_t j = i + 1; j < n; ++j) {
        L(j, i) = A(j, i) / U(i, i);  // Store the factor in L
      }
    }

    return {L, U};
  }

  Decomposition_LDU to_LDU() const {
    Matrix<T> A(*this);
    size_t n = A.nRows();
    Matrix<T> R = A.rowEchelon();  // Row echelon form
    Matrix<T> L(n, n);
    Matrix<T> D(n, n);
    Matrix<T> U = R;  // Initializing U with the row echelon form
    size_t pivot = 1;
    for (size_t i = 0; i < n; ++i) {
      pivot = R(i, i);
      D(i, i) = pivot;
      L(i, i) = 1;

      for (size_t j = i; j < n; j++) {  //normalize U
        U(i, j) /= pivot;
      }

      // Calculate the normalized row echelon form and populate L and U simultaneously
      for (size_t j = i + 1; j < n; ++j) {
        L(j, i) = A(j, i) / U(i, i);
      }
    }

    return Decomposition_LDU(L, D, U);
  }

  std::vector<std::pair<size_t, Matrix<T>>> getPivotColumns() {
    // Perform row echelon form operation
    Matrix<T> rowEchelonMatrix = rowEchelon();

    std::vector<std::pair<size_t, Matrix<T>>> pivotColumns;

    // Iterate through each column of the original matrix
    for (size_t i = 0; i < mCols; ++i) {
      // Check if the corresponding column in the echelon matrix has a pivot
      Matrix<T> Pcolumn = rowEchelonMatrix.getColumn(i);
      bool hasPivot = false;
      for (size_t j = 0; j < Pcolumn.nRows(); ++j) {
        if (Pcolumn(j, 0) != 0) {  // Assuming the matrix stores elements in
                                   // [row][column] format
          hasPivot = true;
          break;
        }
      }

      // If pivot exists, add the corresponding column from the original matrix to
      // the result vector
      if (hasPivot) {
        Matrix<T> originalColumn = getColumn(i);
        pivotColumns.push_back(std::make_pair(i, originalColumn));
      }
    }

    return pivotColumns;
  }

  T det() const {
    try {
      if (nRows() != nCols()) {
        throw not_square();
      }

      auto [echelonMatrix, sign] = echelonizeWithSign();
      T det = 1;  // Initialize determinant value

      for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
        if (echelonMatrix(i, i) == 0) {
          det = 0;
        } else {
          det *= echelonMatrix(i, i);
        }
      }

      return det * sign;
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      throw;
    }
  }
  // end region elementary matrix function

  T Rank() const {
    Matrix<T> rowechelonMat = rowEchelon();

    // Define a tolerance for comparing floating-point numbers
    const double tolerance = 1e-10;

    size_t rank = 0;
    for (size_t i = 0; i < nRows(); i++) {
      Matrix<T> cRow = rowechelonMat.getRow(i);
      // Check if the row is a zero row
      bool isZeroRow = true;
      for (size_t j = 0; j < nCols(); ++j) {
        if (std::abs(cRow(0, j)) > tolerance) {
          isZeroRow = false;
          break;
        }
      }
      if (isZeroRow)
        break;
      else
        rank++;
    }

    return rank;
  }

  Matrix sub_matrix(const size_t row, const size_t column) const {
    if (row >= nRows() || column >= nCols()) {
      throw std::out_of_range("Provided indices are out of bounds.");
    }

    Matrix result(nRows() - 1, nCols() - 1);

    size_t result_row = 0, result_col = 0;

    for (size_t i = 0; i < nRows(); ++i) {
      if (i == row)
        continue;  // Skip the row to be excluded

      for (size_t j = 0; j < nCols(); ++j) {
        if (j == column)
          continue;  // Skip the column to be excluded

        result(result_row, result_col++) = operator()(i, j);
      }

      result_col = 0;  // Reset column index for the next row
      ++result_row;
    }

    return result;
  }

  Matrix transpose() const {
    Matrix result(nRows(), nCols());
    for_ij(nRows(), nCols()) result(j, i) = (*this)(i, j);
    return result;
  }

  Matrix<long double> inverse() const {
    if (mRows != mCols) {
      throw not_square();
    }
    size_t col = nCols();
    // Convert the elements of the original matrix to double
    Matrix<long double> doubleMatrix(*this);
    // Perform inverse calculation with double precision
    Matrix<long double> eye = identity(nRows());
    Matrix<long double> augmented = doubleMatrix.augmentedMatrix(eye);
    Matrix<long double> inverseMatrix =
        augmented.rref()[{{col, 2 * col}, Orientation::Column}];
    return inverseMatrix;
  }

  //------------------------------END CLASS MATRIX--------------------------------------//

};  // class Matrix

// Adjust the number_of_digits function to handle integer values correctly

template <typename T>
size_t number_of_digits(T n) {
  std::ostringstream strs;
  strs << n;
  return strs.str().size();
}

template <typename T>
void print_matrix(const Matrix<T>& matrix) {
  size_t n = matrix.nRows();
  size_t m = matrix.nCols();
  std::vector<size_t> max_len_per_column(m, 0);

  // Calculate the maximum length per column
  for (size_t j = 0; j < m; ++j) {
    for (size_t i = 0; i < n; ++i) {
      size_t num_length = number_of_digits(matrix(i, j));
      if (num_length > max_len_per_column[j])
        max_len_per_column[j] = num_length;
    }
  }

  // Print the matrix
  for (size_t i = 0; i < n; ++i) {
    std::cout << "| ";
    for (size_t j = 0; j < m; ++j) {
      std::cout << std::setw(max_len_per_column[j]) << matrix(i, j);
      if (j != m - 1) {
        std::cout << "   ";
      }
    }
    std::cout << " |\n";
  }
  std::cout << std::endl;
}

// template <typename T>
// Solution<T> solve_AX_b(const Matrix<T>& A, const std::vector<T>& b) {
//     Solution<T> solution;

//     try {
//         // Create the augmented matrix [A|B]
//         Matrix<T> augmented = A.augmentedMatrix(b);

//         // Perform row echelonization
//         Matrix<T> echelonMatrix = augmented.rowEchelon();

//         // Check if the system is consistent (has a solution)
//         for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
//             bool allZero = true;
//             for (size_t j = 0; j < echelonMatrix.nCols() - 1; ++j) {
//                 if (echelonMatrix(i, j) != 0) {
//                     allZero = false;
//                     break;
//                 }
//             }
//             if (allZero && echelonMatrix(i, echelonMatrix.nCols() - 1) != 0) {
//                 // Inconsistent system
//                 solution.type = SolutionType::NO_SOLUTION;
//                 return solution;
//             }
//         }

//         // Initialize the solution vector with zeros
//         solution.values.resize(echelonMatrix.nCols() - 1, T{0});

//         // Backsubstitution to find the solution
//         // Start from the last row and work upwards
//         for (int i = echelonMatrix.nRows() - 1; i >= 0; --i) {
//             T val = echelonMatrix(i, echelonMatrix.nCols() - 1);
//             for (size_t j = i + 1; j < echelonMatrix.nCols() - 1; ++j) {
//                 val -= echelonMatrix(i, j) * solution.values[j];
//             }
//             solution.values[i] = val / echelonMatrix(i, i);
//         }

//         // Determine the solution type
//         if (echelonMatrix.nRows() == echelonMatrix.nCols() - 1) {
//             solution.type = SolutionType::EXACT_SOLUTION;
//         } else {
//             solution.type = SolutionType::INFINITE_SOLUTIONS;
//         }
//     } catch (const std::exception& e) {
//         // Catch any exceptions and print an error message
//         std::cerr << "Error: " << e.what() << std::endl;
//         solution.type = SolutionType::NO_SOLUTION;
//     }

//     return solution;
// }
// Define a template function to solve the equation AX = B
// Define a template function to solve the equation AX = B
// template <typename T>
// Solution<T> solve_AX_b(const Matrix<T>& A, const std::vector<T>& b) {
//     Solution<T> solution;

//     try {
//         // Create the augmented matrix [A|B]
//         Matrix<T> augmented = A.augmentedMatrix(b);

//         // Perform row echelonization
//         Matrix<T> echelonMatrix = augmented.rowEchelon();

//         // Determine the number of leading zeros in each row
//         std::vector<size_t> leadingZeros;
//         for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
//             size_t zeros = 0;
//             while (zeros < echelonMatrix.nCols() - 1 && echelonMatrix(i, zeros) == 0) {
//                 ++zeros;
//             }
//             leadingZeros.push_back(zeros);
//         }

//         // Determine the rank of the matrix
//         size_t rank = 0;
//         for (size_t zeros : leadingZeros) {
//             if (zeros < echelonMatrix.nCols() - 1) {
//                 ++rank;
//             }
//         }

//         // Determine the number of unknowns
//         size_t numUnknowns = echelonMatrix.nCols() - 1;

//         // Determine the number of free variables
//         size_t numFreeVariables = numUnknowns - rank;

//         // Initialize the solution vector
//         solution.values.resize(numUnknowns, T{0});

//         // If there are free variables, set their values to NaN
//         if (numFreeVariables > 0) {
//             for (size_t i = 0; i < numFreeVariables; ++i) {
//                 solution.values[i] = std::numeric_limits<T>::quiet_NaN();
//             }
//         }

//         // Determine the solution type
//         if (rank == numUnknowns) {
//             solution.type = SolutionType::EXACT_SOLUTION;
//         } else if (rank < numUnknowns) {
//             solution.type = SolutionType::INFINITE_SOLUTIONS;
//         } else {
//             solution.type = SolutionType::NO_SOLUTION;
//         }
// if (solution.type == SolutionType::INFINITE_SOLUTIONS) {
//     std::stringstream ss;
//     ss << "Infinite solutions. Linear combination:\n";

//     // For each pivot column
//     for (size_t j = 0; j < numUnknowns; ++j) {
//         if (leadingZeros[j] < echelonMatrix.nRows()) {
//             ss << "x" << j + 1 << " = ";
//             bool firstTerm = true;
//             for (size_t i = 0; i < numUnknowns; ++i) {
//                 if (echelonMatrix(i, j) != 0) {
//                     if (!firstTerm) {
//                         ss << " + ";
//                     } else {
//                         firstTerm = false;
//                     }
//                     if (i < numFreeVariables) {
//                         ss << "c" << i + 1;
//                     } else {
//                         ss << -echelonMatrix(i, j);
//                         for (size_t k = 0; k < numFreeVariables; ++k) {
//                             ss << " * x" << k + 1;
//                             if (k != numFreeVariables - 1) {
//                                 ss << " + ";
//                             }
//                         }
//                     }
//                 }
//             }
//             ss << " + " << echelonMatrix(leadingZeros[j], echelonMatrix.nCols() - 1) << "\n";
//         }
//     }

//     std::cout << ss.str();
// }

//     } catch (const std::exception& e) {
//         // Catch any exceptions and print an error message
//         std::cerr << "Error: " << e.what() << std::endl;
//         solution.type = SolutionType::NO_SOLUTION;
//     }

//     return solution;
// }

// // Function to print the solution
// template<typename T>
// void printSolution(const Solution<T>& solution) {
//     switch (solution.type) {
//         case SolutionType::EXACT_SOLUTION:
//             std::cout << "Exact solution:\n";
//             for (size_t i = 0; i < solution.values.size(); ++i) {
//                 std::cout << "x" << i + 1 << " = " << solution.values[i] << "\n";
//             }
//             break;
//         case SolutionType::INFINITE_SOLUTIONS:
//             // Handled in solve_AX_b function
//             break;
//         case SolutionType::NO_SOLUTION:
//             std::cout << "No solution found.\n";
//             break;
//     }
// }

// Function to generate linear combinations for systems with infinite solutions
// template <typename T>
// std::string linearCombination(const Matrix<T>& echelonMatrix) {
//     std::string combination;
//     std::vector<bool> freeVariable(echelonMatrix.nCols() - 1, true);

//     // Mark the pivot columns
//     for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
//         for (size_t j = 0; j < echelonMatrix.nCols() - 1; ++j) {
//             if (echelonMatrix(i, j) != 0) {
//                 freeVariable[j] = false;
//                 break;
//             }
//         }
//     }

//     for (size_t i = 0; i < freeVariable.size(); ++i) {
//         if (freeVariable[i]) {
//             std::string variable = "x" + std::to_string(i + 1) + " = ";
//             bool firstTerm = true;
//             for (size_t j = 0; j < echelonMatrix.nCols() - 1; ++j) {
//                 if (j != i && echelonMatrix(j, i) != 0) {
//                     if (!firstTerm) {
//                         variable += " + ";
//                     }
//                     if (echelonMatrix(j, i) == 1) {
//                         variable += "c" + std::to_string(j + 1);
//                     } else if (echelonMatrix(j, i) == -1) {
//                         variable += "-c" + std::to_string(j + 1);
//                     } else {
//                         variable += std::to_string(echelonMatrix(j, i)) + "c" + std::to_string(j + 1);
//                     }
//                     firstTerm = false;
//                 }
//             }
//             combination += variable + "\n";
//         }
//     }

//     return combination;
// }

// template <typename T>
// Solution<T> solve_AX_b(const Matrix<T>& A, const std::vector<T>& b) {
//     Solution<T> solution;

//     // Get the row echelon form
//     Matrix<T> echelonMatrix = A.rowEchelon();
//     std::cout << "Matrix A and b augmented row echelon:\n";
//     print_matrix(echelonMatrix);

//     // Check if the system is consistent (has a solution)
//     for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
//         bool allZero = true;
//         for (size_t j = 0; j < echelonMatrix.nCols() - 1; ++j) {
//             if (echelonMatrix(i, j) != 0) {
//                 allZero = false;
//                 break;
//             }
//         }
//         if (allZero && echelonMatrix(i, echelonMatrix.nCols() - 1) != 0) {
//             solution.type = SolutionType::NO_SOLUTION;
//             return solution;
//         }
//     }

//     // Check if the system has unique solution
//     if (echelonMatrix.nRows() == echelonMatrix.nCols() - 1) {
//         solution.type = SolutionType::EXACT_SOLUTION;
//         solution.values = back_substitution(echelonMatrix);
//         return solution;
//     }

//     // Check if the system has infinite solutions
//     if (echelonMatrix.nRows() < echelonMatrix.nCols() - 1) {
//         solution.type = SolutionType::INFINITE_SOLUTIONS;
//         solution.values.clear(); // Clear existing values
//         solution.values.push_back("Infinite solutions. Linear combination:");
//         solution.values.push_back(linearCombination(echelonMatrix));
//         return solution;
//     }

//     // If none of the above conditions are met, return no solution
//     solution.type = SolutionType::NO_SOLUTION;
//     return solution;
// }

// Function to generate linear combinattemplate <typename T>

// Function to generate linear combinations for systems with infinite solutions

// Function to solve the system AX = b
// template <typename T>
// Solution<T> solve_AX_b(const Matrix<T>& A, const std::vector<T>& b) {
//     Solution<T> solution;

//     // Get the row echelon form
//     Matrix<T> echelonMatrix = A.augmentedMatrix(b).rowEchelon();
//     std::cout << "Matrix A and b augmented row echelon:\n";
//     print_matrix(echelonMatrix);

//     // Calculate the rank of the augmented matrix
//     size_t rank = echelonMatrix.Rank();

//     // Determine the number of variables
//     size_t numVariables = echelonMatrix.nCols() - 1;

//     // Check if the system is consistent (has a solution)
//     if (!isConsistent(echelonMatrix)) {
//         solution.type = SolutionType::NO_SOLUTION;
//         return solution;
//     }

//     // Check if the system has unique solution
//     if (rank == numVariables) {
//         solution.type = SolutionType::EXACT_SOLUTION;
//         solution.values = backSubstitution(echelonMatrix);
//         return solution;
//     }

//     // Check if the system has infinite solutions
//     if (rank < numVariables) {
//         solution.type = SolutionType::INFINITE_SOLUTIONS;
//         solution.linearComb = linearCombination(echelonMatrix);
//         return solution;
//     }

//     return solution; // Shouldn't reach here, but return empty solution as fallback
// }

enum class SolutionType {
  EXACT_SOLUTION_XP,         // Unique solution exists
  INFINITE_SOLUTIONS_XP_XS,  // Infinite solutions exist
  NO_SOLUTION,
  ERROR  // No solution exists
};

template <typename T>
struct Solution {
  SolutionType type;

  std::vector<T> values;
  std::string linearComb;  // For infinite solutions
};

template <typename T>
bool isConsistent(const Matrix<T>& echelonMatrix) {
  for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
    bool allZero = true;
    for (size_t j = 0; j < echelonMatrix.nCols() - 1; ++j) {
      if (echelonMatrix(i, j) != 0) {
        allZero = false;
        break;
      }
    }
    if (allZero && echelonMatrix(i, echelonMatrix.nCols() - 1) != 0) {
      return false;  // Inconsistent equations
    }
  }
  return true;  // Consistent equations
}

// Function to check if the system has a unique solution
template <typename T>
bool hasUniqueSolution(const Matrix<T>& echelonMatrix) {
  return echelonMatrix.nRows() == echelonMatrix.nCols() - 1;
}

// Function to perform back substitution
template <typename T>
std::vector<T> backSubstitution(const Matrix<T>& echelonMatrix) {
  std::vector<T> solution(echelonMatrix.nRows(), 0);

  for (int i = echelonMatrix.nRows() - 1; i >= 0; --i) {
    T val = echelonMatrix(i, echelonMatrix.nCols() - 1);
    for (size_t j = i + 1; j < echelonMatrix.nCols() - 1; ++j) {
      val -= echelonMatrix(i, j) * solution[j];
    }
    solution[i] = val / echelonMatrix(i, i);
  }

  return solution;
}

template <typename T>
bool hasZeroPivotColumn(const Matrix<T>& echelonMatrix) {
  for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
    if (echelonMatrix(i, i) == 0) {
      return true;  // At least one diagonal element (pivot) is zero
    }
  }
  return false;  // All diagonal elements (pivots) are non-zero
}

// template <typename T>
// bool isZeroRow(const Matrix<T>& matrix, size_t rowIdx) {
//   size_t numCols = matrix.nCols();
//   for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
//     if (matrix(rowIdx, colIdx) != 0) {
//       return false;
//     }
//   }
//   return true;
// }
template <typename T>
void findPivotAndFreeColumns(const Matrix<T>& rrefMatrix,
                             std::vector<size_t>& pivotColumns,
                             std::vector<size_t>& freeColumns) {
  size_t numRows = rrefMatrix.nRows();
  size_t numCols = rrefMatrix.nCols() - 1;  // Exclude the last column (b)

  for (size_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
    bool foundPivot = false;
    for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
      if (rrefMatrix(rowIdx, colIdx) != 0) {
        pivotColumns.push_back(colIdx);
        foundPivot = true;
        break;
      }
    }
    if (!foundPivot) {
      // If no pivot is found in this row, all columns in this row are free columns
      for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
        freeColumns.push_back(colIdx);
      }
      // No need to check further rows
      break;
    }
  }

  // Remove pivot columns from freeColumns
  std::sort(pivotColumns.begin(), pivotColumns.end());
  std::sort(freeColumns.begin(), freeColumns.end());
  std::vector<size_t> intersection;
  std::set_intersection(pivotColumns.begin(), pivotColumns.end(),
                        freeColumns.begin(), freeColumns.end(),
                        std::back_inserter(intersection));
  for (auto col : intersection) {
    freeColumns.erase(std::remove(freeColumns.begin(), freeColumns.end(), col),
                      freeColumns.end());
  }
}
template <typename T>
void calculateXP(const Matrix<T>& rrefMatrix,
                 const std::vector<size_t>& pivotColumns, size_t numUnknowns,
                 std::vector<T>& xp) {
  for (size_t i = 0; i < pivotColumns.size(); ++i) {
    size_t pivotCol = pivotColumns[i];
    for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
      if (rrefMatrix(rowIdx, pivotCol) != 0) {
        xp[pivotCol] = rrefMatrix(rowIdx, numUnknowns);
        break;
      }
    }
  }
}
template <typename T>
std::vector<std::string> calculateXS(const Matrix<T>& rrefMatrix,
                                     const std::vector<size_t>& pivotColumns,
                                     const std::vector<size_t>& freeColumns,
                                     size_t numCols) {
  std::vector<std::string> xsExpressions;
  size_t numXs = freeColumns.size();

  for (size_t xsIdx = 0; xsIdx < numXs; ++xsIdx) {
    std::vector<T> xs(numCols, 0);
    size_t freeCol = freeColumns[xsIdx];
    xs[freeCol] = 1;

    for (size_t i = 0; i < pivotColumns.size(); ++i) {
      size_t pivotCol = pivotColumns[i];
      bool foundNonZero = false;
      T coefficient = 0;  // Initialize coefficient to 0
      for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
        if (rrefMatrix(rowIdx, pivotCol) != 0) {
          if (!foundNonZero) {
            coefficient =
                -rrefMatrix(rowIdx, freeCol) / rrefMatrix(rowIdx, pivotCol);
            foundNonZero = true;
          } else {
            coefficient *= rrefMatrix(rowIdx, pivotCol);
            coefficient -= rrefMatrix(rowIdx, freeCol);
          }
        }
      }
      xs[pivotCol] = coefficient;
    }

    std::ostringstream xsStream;
    xsStream << std::fixed << std::setprecision(2);
    xsStream << "Xs" << xsExpressions.size() + 1 << ": ";
    for (size_t i = 0; i < numCols; ++i) {
      xsStream << "X" << i + 1 << " = " << xs[i];
      if (i < numCols - 1)
        xsStream << "; ";
    }
    xsExpressions.push_back(xsStream.str());
  }

  return xsExpressions;
}

std::string LinearCombination(const std::vector<std::string>& xsExpressions) {
  std::string linearComb;
  for (size_t i = 0; i < xsExpressions.size(); ++i) {
    linearComb += "C" + std::to_string(i + 1) + " * (" + xsExpressions[i] + ")";
    if (i < xsExpressions.size() - 1)
      linearComb += " + ";
  }
  return linearComb;
}

// Your solver_AX_b function
template <typename T>
Solution<T> solver_AX_b(const Matrix<T>& A, const std::vector<T>& b) {
  Matrix<T> augmented = A.augmentedMatrix(b);

  // Print row echelon form of augmented matrix for debugging exact solutions
  std::cout << "Row echelon form of augmented matrix:" << std::endl;
  Matrix<T> echelonMatrix = augmented.rowEchelon();

  print_matrix(echelonMatrix);  // Assuming print_matrix is defined

  // Calculate the rank of the augmented matrix
  size_t rank = echelonMatrix.Rank();
  size_t numUnknowns = echelonMatrix.nCols() - 1;

  if (rank == echelonMatrix.nCols() - 1) {
    // If rank equals the number of unknowns, it has an exact solution
    std::vector<T> xp = backSubstitution(echelonMatrix);
    return {SolutionType::EXACT_SOLUTION_XP, xp, ""};
  } else {
    // Check if there is at least one row where the last column is zero
    bool noSolution = true;
    for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
      if (echelonMatrix(i, echelonMatrix.nCols() - 1) == 0) {
        noSolution = false;
        break;
      }
    }

    if (noSolution) {
      return {SolutionType::NO_SOLUTION, {}, "No solution exists."};
    } else {
      // Determine special solutions for infinite solutions
      Matrix<T> rrefMatrix = augmented.rref();
      size_t numRows = rrefMatrix.nRows();
      size_t numCols = rrefMatrix.nCols() - 1;
      std::vector<size_t> pivotColumns;
      std::vector<size_t> freeColumns;

      findPivotAndFreeColumns(rrefMatrix, pivotColumns, freeColumns);

      std::vector<T> xp(numCols, 0);
      calculateXP(rrefMatrix, pivotColumns, numUnknowns, xp);

      // Calculate special solutions XS
      std::vector<std::string> xsExpressions =
          calculateXS(rrefMatrix, pivotColumns, freeColumns, numCols);

      // Generate linear combination string
      std::string linearComb = LinearCombination(xsExpressions);

      return {SolutionType::INFINITE_SOLUTIONS_XP_XS, xp, linearComb};
    }
  }
}

template <typename T>
Solution<T> solver_AX_b33(const Matrix<T>& A, const std::vector<T>& b) {
  Matrix<T> augmented = A.augmentedMatrix(b);

  // Print row echelon form of augmented matrix for debugging exact solutions
  std::cout << "Row echelon form of augmented matrix:" << std::endl;
  Matrix<T> echelonMatrix = augmented.rowEchelon();
  print_matrix(echelonMatrix);

  // Calculate the rank of the augmented matrix

  size_t rank = echelonMatrix.Rank();
  size_t numUnknowns = echelonMatrix.nCols() - 1;

  if (rank == echelonMatrix.nCols() - 1) {
    // If rank equals the number of unknowns, it has an exact solution
    std::vector<T> xp = backSubstitution(echelonMatrix);
    return {SolutionType::EXACT_SOLUTION_XP, xp, ""};
  } else {
    // Check if there is at least one row where the last column is zero
    bool noSolution = true;
    for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
      if (echelonMatrix(i, echelonMatrix.nCols() - 1) == 0) {
        noSolution = false;
        break;
      }
    }

    if (noSolution) {
      return {SolutionType::NO_SOLUTION, {}, "No solution exists."};
    } else {
      // Determine special solutions for infinite solutions
      Matrix<T> rrefMatrix = augmented.rref();
      std::vector<size_t> pivotColumns;
      std::vector<size_t> freeColumns;
      size_t numRows = rrefMatrix.nRows();
      size_t numCols = rrefMatrix.nCols() - 1;
      findPivotAndFreeColumns(rrefMatrix, pivotColumns, freeColumns);

      std::vector<T> xp(numCols, 0);

      calculateXP(rrefMatrix, pivotColumns, numUnknowns, xp);
      std::vector<std::string> xsExpressions =
          calculateXS(rrefMatrix, pivotColumns, freeColumns, numCols);

      std::string linearComb;
      for (size_t i = 0; i < xsExpressions.size(); ++i) {
        linearComb +=
            "C" + std::to_string(i + 1) + " * (" + xsExpressions[i] + ")";
        if (i < xsExpressions.size() - 1)
          linearComb += " + ";
      }

      return {SolutionType::INFINITE_SOLUTIONS_XP_XS, xp, linearComb};
    }
  }
}

template <typename T>
Solution<T> solver_AX_b22(const Matrix<T>& A, const std::vector<T>& b) {
  Matrix<T> augmented = A.augmentedMatrix(b);

  // Print row echelon form of augmented matrix for debugging exact solutions
  std::cout << "Row echelon form of augmented matrix:" << std::endl;
  Matrix<T> echelonMatrix = augmented.rowEchelon();
  print_matrix(echelonMatrix);

  // Calculate the rank of the augmented matrix
  size_t rank = echelonMatrix.Rank();

  if (rank == echelonMatrix.nCols() - 1) {
    // If rank equals the number of unknowns, it has an exact solution
    std::vector<T> xp = backSubstitution(echelonMatrix);
    return {SolutionType::EXACT_SOLUTION_XP, xp, ""};
  } else {
    // Check if there is at least one row where the last column is zero
    bool noSolution = true;
    for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
      if (echelonMatrix(i, echelonMatrix.nCols() - 1) == 0) {
        noSolution = false;
        break;
      }
    }

    if (noSolution) {
      return {SolutionType::NO_SOLUTION, {}, "No solution exists."};
    } else {
      // Determine special solutions for infinite solutions
      Matrix<T> rrefMatrix = augmented.rref();
      std::vector<size_t> pivotColumns;
      std::vector<size_t> freeColumns;
      size_t numRows = rrefMatrix.nRows();
      size_t numCols = rrefMatrix.nCols() - 1;
      findPivotAndFreeColumns(rrefMatrix, pivotColumns, freeColumns);

      std::vector<T> xp(numCols, 0);

      for (size_t i = 0; i < pivotColumns.size(); ++i) {
        size_t pivotCol = pivotColumns[i];
        for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
          if (rrefMatrix(rowIdx, pivotCol) != 0) {
            xp[pivotCol] = rrefMatrix(rowIdx, numCols);
            break;
          }
        }
      }

      std::vector<std::string> xsExpressions;
      size_t numFreeVars = freeColumns.size();
      size_t numXs = numFreeVars;
      for (size_t xsIdx = 0; xsIdx < numXs; ++xsIdx) {
        std::vector<T> xs(numCols, 0);
        size_t freeCol = freeColumns[xsIdx];
        xs[freeCol] = 1;

        for (size_t i = 0; i < pivotColumns.size(); ++i) {
          size_t pivotCol = pivotColumns[i];
          bool foundNonZero = false;
          T coefficient;
          for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
            if (rrefMatrix(rowIdx, pivotCol) != 0) {
              if (!foundNonZero) {
                coefficient =
                    -rrefMatrix(rowIdx, freeCol) / rrefMatrix(rowIdx, pivotCol);
                foundNonZero = true;
              } else {
                coefficient *= rrefMatrix(rowIdx, pivotCol);
                coefficient -= rrefMatrix(rowIdx, freeCol);
              }
            }
          }
          xs[pivotCol] = coefficient;
        }

        std::ostringstream xsStream;
        xsStream << std::fixed << std::setprecision(2);
        xsStream << "Xs" << xsExpressions.size() + 1 << ": ";
        for (size_t i = 0; i < numCols; ++i) {
          xsStream << "X" << i + 1 << " = " << xs[i];
          if (i < numCols - 1)
            xsStream << "; ";
        }
        xsExpressions.push_back(xsStream.str());
      }

      std::string linearComb;
      for (size_t i = 0; i < xsExpressions.size(); ++i) {
        linearComb +=
            "C" + std::to_string(i + 1) + " * (" + xsExpressions[i] + ")";
        if (i < xsExpressions.size() - 1)
          linearComb += " + ";
      }

      return {SolutionType::INFINITE_SOLUTIONS_XP_XS, xp, linearComb};
    }
  }
}

template <typename T>
Solution<T> solver_AX_b11(const Matrix<T>& A, const std::vector<T>& b) {
  Matrix<T> augmented = A.augmentedMatrix(b);

  // Print row echelon form of augmented matrix for debugging exact solutions
  std::cout << "Row echelon form of augmented matrix:" << std::endl;
  Matrix<T> echelonMatrix = augmented.rowEchelon();
  print_matrix(echelonMatrix);

  // Calculate the rank of the augmented matrix
  size_t rank = echelonMatrix.Rank();

  if (rank == echelonMatrix.nCols() - 1) {
    // If rank equals the number of unknowns, it has an exact solution
    std::vector<T> xp = backSubstitution(echelonMatrix);
    return {SolutionType::EXACT_SOLUTION_XP, xp, ""};
  } else {
    // Check if there is at least one row where the last column is zero
    bool noSolution = true;
    for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
      if (echelonMatrix(i, echelonMatrix.nCols() - 1) == 0) {
        noSolution = false;
        break;
      }
    }

    if (noSolution) {
      return {SolutionType::NO_SOLUTION, {}, "No solution exists."};
    } else {
      // Determine special solutions for infinite solutions
      Matrix<T> rrefMatrix = augmented.rref();
      std::vector<size_t> pivotColumns;
      std::vector<size_t> freeColumns;
      size_t numRows = rrefMatrix.nRows();
      size_t numCols = rrefMatrix.nCols() - 1;  // Exclude the last column (b)

      // Iterate through rows to find pivot columns
      for (size_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
        bool foundPivot = false;
        for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
          if (rrefMatrix(rowIdx, colIdx) != 0) {
            pivotColumns.push_back(colIdx);
            foundPivot = true;
            break;
          }
        }
        if (!foundPivot) {
          // If no pivot is found in this row, all columns in this row are free columns
          for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
            freeColumns.push_back(colIdx);
          }
          // No need to check further rows
          break;
        }
      }

      // Remove pivot columns from freeColumns
      std::sort(pivotColumns.begin(), pivotColumns.end());
      std::sort(freeColumns.begin(), freeColumns.end());
      std::vector<size_t> intersection;
      std::set_intersection(pivotColumns.begin(), pivotColumns.end(),
                            freeColumns.begin(), freeColumns.end(),
                            std::back_inserter(intersection));
      for (auto col : intersection) {
        freeColumns.erase(
            std::remove(freeColumns.begin(), freeColumns.end(), col),
            freeColumns.end());
      }

      std::vector<T> xp(numCols, 0);

      for (size_t i = 0; i < pivotColumns.size(); ++i) {
        size_t pivotCol = pivotColumns[i];
        for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
          if (rrefMatrix(rowIdx, pivotCol) != 0) {
            xp[pivotCol] = rrefMatrix(rowIdx, numCols);
            break;
          }
        }
      }

      std::vector<std::string> xsExpressions;
      size_t numFreeVars = freeColumns.size();
      size_t numXs = numFreeVars;
      for (size_t xsIdx = 0; xsIdx < numXs; ++xsIdx) {
        std::vector<T> xs(numCols, 0);
        size_t freeCol = freeColumns[xsIdx];
        xs[freeCol] = 1;

        for (size_t i = 0; i < pivotColumns.size(); ++i) {
          size_t pivotCol = pivotColumns[i];
          bool foundNonZero = false;
          T coefficient;
          for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
            if (rrefMatrix(rowIdx, pivotCol) != 0) {
              if (!foundNonZero) {
                coefficient =
                    -rrefMatrix(rowIdx, freeCol) / rrefMatrix(rowIdx, pivotCol);
                foundNonZero = true;
              } else {
                coefficient *= rrefMatrix(rowIdx, pivotCol);
                coefficient -= rrefMatrix(rowIdx, freeCol);
              }
            }
          }
          xs[pivotCol] = coefficient;
        }

        std::ostringstream xsStream;
        xsStream << std::fixed << std::setprecision(2);
        xsStream << "Xs" << xsExpressions.size() + 1 << ": ";
        for (size_t i = 0; i < numCols; ++i) {
          xsStream << "X" << i + 1 << " = " << xs[i];
          if (i < numCols - 1)
            xsStream << "; ";
        }
        xsExpressions.push_back(xsStream.str());
      }

      std::string linearComb;
      for (size_t i = 0; i < xsExpressions.size(); ++i) {
        linearComb +=
            "C" + std::to_string(i + 1) + " * (" + xsExpressions[i] + ")";
        if (i < xsExpressions.size() - 1)
          linearComb += " + ";
      }

      return {SolutionType::INFINITE_SOLUTIONS_XP_XS, xp, linearComb};
    }
  }
}

template <typename T>
Solution<T> solver_AX_b4(const Matrix<T>& A, const std::vector<T>& b) {
  Matrix<T> augmented = A.augmentedMatrix(b);

  // Print row echelon form of augmented matrix for debugging exact solutions
  std::cout << "Row echelon form of augmented matrix:" << std::endl;
  Matrix<T> echelonMatrix = augmented.rowEchelon();
  print_matrix(echelonMatrix);

  // Check if all pivots are non-zero to determine exact or infinite solution
  bool allNonZeroPivots = true;
  for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
    if (echelonMatrix(i, i) == 0) {
      allNonZeroPivots = false;
      break;
    }
  }

  if (allNonZeroPivots) {
    // If all pivots are non-zero, calculate exact solution (Xp) using back substitution
    std::vector<T> xp = backSubstitution(echelonMatrix);
    return {SolutionType::EXACT_SOLUTION_XP, xp, ""};
  } else {
    // If one or more pivots are zero, indicating infinite solutions, Xp should have pivot variables
    Matrix<T> rrefMatrix = augmented.rref();

    std::cout << "Reduced row echelon form of augmented matrix:" << std::endl;
    print_matrix(rrefMatrix);

    std::vector<size_t> pivotColumns;
    std::vector<size_t> freeColumns;
    size_t numRows = rrefMatrix.nRows();
    size_t numCols = rrefMatrix.nCols() - 1;  // Exclude the last column (b)

    // Iterate through rows to find pivot columns
    for (size_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
      bool foundPivot = false;
      for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
        if (rrefMatrix(rowIdx, colIdx) != 0) {
          pivotColumns.push_back(colIdx);
          foundPivot = true;
          break;
        }
      }
      if (!foundPivot) {
        // If no pivot is found in this row, all columns in this row are free columns
        for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
          freeColumns.push_back(colIdx);
        }
        // No need to check further rows
        break;
      }
    }

    // If the last column does not contain a pivot column, the system is inconsistent
    if (pivotColumns.empty() || pivotColumns.back() != numCols - 1) {
      return {SolutionType::NO_SOLUTION, {}, "No solution exists."};
    }

    // Remove pivot columns from freeColumns
    std::sort(pivotColumns.begin(), pivotColumns.end());
    std::sort(freeColumns.begin(), freeColumns.end());
    std::vector<size_t> intersection;
    std::set_intersection(pivotColumns.begin(), pivotColumns.end(),
                          freeColumns.begin(), freeColumns.end(),
                          std::back_inserter(intersection));
    for (auto col : intersection) {
      freeColumns.erase(
          std::remove(freeColumns.begin(), freeColumns.end(), col),
          freeColumns.end());
    }

    std::cout << "Pivot columns: ";
    for (size_t i = 0; i < pivotColumns.size(); ++i) {
      std::cout << pivotColumns[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Free columns: ";
    for (size_t i = 0; i < freeColumns.size(); ++i) {
      std::cout << freeColumns[i] << " ";
    }
    std::cout << std::endl;

    std::vector<T> xp(numCols, 0);

    for (size_t i = 0; i < pivotColumns.size(); ++i) {
      size_t pivotCol = pivotColumns[i];
      for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
        if (rrefMatrix(rowIdx, pivotCol) != 0) {
          xp[pivotCol] = rrefMatrix(rowIdx, numCols);
          break;
        }
      }
    }

    std::vector<std::string> xsExpressions;
    size_t numFreeVars = freeColumns.size();
    size_t numXs = numFreeVars;
    for (size_t xsIdx = 0; xsIdx < numXs; ++xsIdx) {
      std::vector<T> xs(numCols, 0);
      size_t freeCol = freeColumns[xsIdx];
      xs[freeCol] = 1;

      // Adjust xs values for special solutions
      for (size_t i = 0; i < pivotColumns.size(); ++i) {
        size_t pivotCol = pivotColumns[i];
        bool foundNonZero = false;
        T coefficient;
        for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
          if (rrefMatrix(rowIdx, pivotCol) != 0) {
            if (!foundNonZero) {
              coefficient =
                  -rrefMatrix(rowIdx, freeCol) / rrefMatrix(rowIdx, pivotCol);
              foundNonZero = true;
            } else {
              coefficient *=
                  -rrefMatrix(rowIdx, freeCol) / rrefMatrix(rowIdx, pivotCol);
            }
          }
        }
        xs[pivotCol] = coefficient;
      }

      std::ostringstream xsStream;
      xsStream << std::fixed << std::setprecision(2);
      xsStream << "Xs" << xsExpressions.size() + 1 << ": ";
      for (size_t i = 0; i < numCols; ++i) {
        xsStream << "X" << i + 1 << " = " << xs[i];
        if (i < numCols - 1)
          xsStream << "; ";
      }
      xsExpressions.push_back(xsStream.str());
    }

    std::string linearComb;
    for (size_t i = 0; i < xsExpressions.size(); ++i) {
      linearComb +=
          "C" + std::to_string(i + 1) + " * (" + xsExpressions[i] + ")";
      if (i < xsExpressions.size() - 1)
        linearComb += " + ";
    }

    return {SolutionType::INFINITE_SOLUTIONS_XP_XS, xp, linearComb};
  }
}

template <typename T>
Solution<T> solver_AX_b3(const Matrix<T>& A, const std::vector<T>& b) {
  Matrix<T> augmented = A.augmentedMatrix(b);

  // Print row echelon form of augmented matrix for debugging exact solutions
  std::cout << "Row echelon form of augmented matrix:" << std::endl;
  Matrix<T> echelonMatrix = augmented.rowEchelon();
  print_matrix(echelonMatrix);

  // Check if all pivots are non-zero to determine exact or infinite solution
  bool allNonZeroPivots = true;
  for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
    if (echelonMatrix(i, i) == 0) {
      allNonZeroPivots = false;
      break;
    }
  }

  if (allNonZeroPivots) {
    // If all pivots are non-zero, calculate exact solution (Xp) using back substitution
    std::vector<T> xp = backSubstitution(echelonMatrix);
    return {SolutionType::EXACT_SOLUTION_XP, xp, ""};
  } else {
    // If one or more pivots are zero, indicating infinite solutions, Xp should have pivot variables
    Matrix<T> rrefMatrix = augmented.rref();

    std::cout << "Reduced row echelon form of augmented matrix:" << std::endl;
    print_matrix(rrefMatrix);

    std::vector<size_t> pivotColumns;
    std::vector<size_t> freeColumns;
    size_t numRows = rrefMatrix.nRows();
    size_t numCols = rrefMatrix.nCols() - 1;  // Exclude the last column (b)

    // Iterate through rows to find pivot columns
    for (size_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
      bool foundPivot = false;
      for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
        if (rrefMatrix(rowIdx, colIdx) != 0) {
          pivotColumns.push_back(colIdx);
          foundPivot = true;
          break;
        }
      }
      if (!foundPivot) {
        // If no pivot is found in this row, all columns in this row are free columns
        for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
          freeColumns.push_back(colIdx);
        }
        // No need to check further rows
        break;
      }
    }

    // Remove pivot columns from freeColumns
    std::sort(pivotColumns.begin(), pivotColumns.end());
    std::sort(freeColumns.begin(), freeColumns.end());
    std::vector<size_t> intersection;
    std::set_intersection(pivotColumns.begin(), pivotColumns.end(),
                          freeColumns.begin(), freeColumns.end(),
                          std::back_inserter(intersection));
    for (auto col : intersection) {
      freeColumns.erase(
          std::remove(freeColumns.begin(), freeColumns.end(), col),
          freeColumns.end());
    }

    std::cout << "Pivot columns: ";
    for (size_t i = 0; i < pivotColumns.size(); ++i) {
      std::cout << pivotColumns[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Free columns: ";
    for (size_t i = 0; i < freeColumns.size(); ++i) {
      std::cout << freeColumns[i] << " ";
    }
    std::cout << std::endl;

    std::vector<T> xp(numCols, 0);

    for (size_t i = 0; i < pivotColumns.size(); ++i) {
      size_t pivotCol = pivotColumns[i];
      for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
        if (rrefMatrix(rowIdx, pivotCol) != 0) {
          xp[pivotCol] = rrefMatrix(rowIdx, numCols);
          break;
        }
      }
    }

    std::vector<std::string> xsExpressions;
    size_t numFreeVars = freeColumns.size();
    size_t numXs = numFreeVars;
    for (size_t xsIdx = 0; xsIdx < numXs; ++xsIdx) {
      std::vector<T> xs(numCols, 0);
      size_t freeCol = freeColumns[xsIdx];
      xs[freeCol] = 1;

      // Adjust xs values for special solutions
      for (size_t i = 0; i < pivotColumns.size(); ++i) {
        size_t pivotCol = pivotColumns[i];
        bool foundNonZero = false;
        T coefficient;
        for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
          if (rrefMatrix(rowIdx, pivotCol) != 0) {
            if (!foundNonZero) {
              coefficient =
                  -rrefMatrix(rowIdx, freeCol) / rrefMatrix(rowIdx, pivotCol);
              foundNonZero = true;
            } else {
              coefficient *=
                  -rrefMatrix(rowIdx, freeCol) / rrefMatrix(rowIdx, pivotCol);
            }
          }
        }
        xs[pivotCol] = coefficient;
      }

      std::ostringstream xsStream;
      xsStream << std::fixed << std::setprecision(2);
      xsStream << "Xs" << xsExpressions.size() + 1 << ": ";
      for (size_t i = 0; i < numCols; ++i) {
        xsStream << "X" << i + 1 << " = " << xs[i];
        if (i < numCols - 1)
          xsStream << "; ";
      }
      xsExpressions.push_back(xsStream.str());
    }

    std::string linearComb;
    for (size_t i = 0; i < xsExpressions.size(); ++i) {
      linearComb +=
          "C" + std::to_string(i + 1) + " * (" + xsExpressions[i] + ")";
      if (i < xsExpressions.size() - 1)
        linearComb += " + ";
    }

    return {SolutionType::INFINITE_SOLUTIONS_XP_XS, xp, linearComb};
  }
}

//correct one
template <typename T>
Solution<T> solver_AX_b2(const Matrix<T>& A, const std::vector<T>& b) {
  Matrix<T> augmented = A.augmentedMatrix(b);

  // Print row echelon form of augmented matrix for debugging exact solutions
  std::cout << "Row echelon form of augmented matrix:" << std::endl;
  Matrix<T> echelonMatrix = augmented.rowEchelon();
  print_matrix(echelonMatrix);

  // Check if all pivots are non-zero to determine exact or infinite solution
  bool allNonZeroPivots = true;
  for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
    if (echelonMatrix(i, i) == 0) {
      allNonZeroPivots = false;
      break;
    }
  }

  if (allNonZeroPivots) {
    // If all pivots are non-zero, calculate exact solution (Xp) using back substitution
    std::vector<T> xp = backSubstitution(echelonMatrix);
    return {SolutionType::EXACT_SOLUTION_XP, xp, ""};
  } else {
    // If one or more pivots are zero, indicating infinite solutions, Xp should have pivot variables
    Matrix<T> rrefMatrix = augmented.rref();

    std::cout << "Reduced row echelon form of augmented matrix:" << std::endl;
    print_matrix(rrefMatrix);

    std::vector<size_t> pivotColumns;
    std::vector<size_t> freeColumns;
    size_t numRows = rrefMatrix.nRows();
    size_t numCols = rrefMatrix.nCols() - 1;  // Exclude the last column (b)

    // Iterate through rows to find pivot columns
    for (size_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
      bool foundPivot = false;
      for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
        if (rrefMatrix(rowIdx, colIdx) != 0) {
          pivotColumns.push_back(colIdx);
          foundPivot = true;
          break;
        }
      }
      if (!foundPivot) {
        // If no pivot is found in this row, all columns in this row are free columns
        for (size_t colIdx = 0; colIdx < numCols; ++colIdx) {
          freeColumns.push_back(colIdx);
        }
        // No need to check further rows
        break;
      }
    }

    // Remove pivot columns from freeColumns
    std::sort(pivotColumns.begin(), pivotColumns.end());
    std::sort(freeColumns.begin(), freeColumns.end());
    std::vector<size_t> intersection;
    std::set_intersection(pivotColumns.begin(), pivotColumns.end(),
                          freeColumns.begin(), freeColumns.end(),
                          std::back_inserter(intersection));
    for (auto col : intersection) {
      freeColumns.erase(
          std::remove(freeColumns.begin(), freeColumns.end(), col),
          freeColumns.end());
    }

    std::cout << "Pivot columns: ";
    for (size_t i = 0; i < pivotColumns.size(); ++i) {
      std::cout << pivotColumns[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Free columns: ";
    for (size_t i = 0; i < freeColumns.size(); ++i) {
      std::cout << freeColumns[i] << " ";
    }
    std::cout << std::endl;

    std::vector<T> xp(numCols, 0);

    for (size_t i = 0; i < pivotColumns.size(); ++i) {
      size_t pivotCol = pivotColumns[i];
      for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
        if (rrefMatrix(rowIdx, pivotCol) != 0) {
          xp[pivotCol] = rrefMatrix(rowIdx, numCols);
          break;
        }
      }
    }

    std::vector<std::string> xsExpressions;
    size_t numFreeVars = freeColumns.size();
    size_t numXs = numFreeVars;
    for (size_t xsIdx = 0; xsIdx < numXs; ++xsIdx) {
      std::vector<T> xs(numCols, 0);
      size_t freeCol = freeColumns[xsIdx];
      xs[freeCol] = 1;

      // Adjust xs values for special solutions
      for (size_t i = 0; i < pivotColumns.size(); ++i) {
        size_t pivotCol = pivotColumns[i];
        bool foundNonZero = false;
        T coefficient;
        for (size_t rowIdx = 0; rowIdx < rrefMatrix.nRows(); ++rowIdx) {
          if (rrefMatrix(rowIdx, pivotCol) != 0) {
            if (!foundNonZero) {
              coefficient =
                  -rrefMatrix(rowIdx, freeCol) / rrefMatrix(rowIdx, pivotCol);
              foundNonZero = true;
            } else {
              coefficient *= rrefMatrix(rowIdx, pivotCol);
              coefficient -= rrefMatrix(rowIdx, freeCol);
            }
          }
        }
        xs[pivotCol] = coefficient;
      }

      std::ostringstream xsStream;
      xsStream << std::fixed << std::setprecision(2);
      xsStream << "Xs" << xsExpressions.size() + 1 << ": ";
      for (size_t i = 0; i < numCols; ++i) {
        xsStream << "X" << i + 1 << " = " << xs[i];
        if (i < numCols - 1)
          xsStream << "; ";
      }
      xsExpressions.push_back(xsStream.str());
    }

    std::string linearComb;
    for (size_t i = 0; i < xsExpressions.size(); ++i) {
      linearComb +=
          "C" + std::to_string(i + 1) + " * (" + xsExpressions[i] + ")";
      if (i < xsExpressions.size() - 1)
        linearComb += " + ";
    }

    return {SolutionType::INFINITE_SOLUTIONS_XP_XS, xp, linearComb};
  }
}

// template <typename T>
// Solution<T> solver_AX_b(const Matrix<T>& A, const std::vector<T>& b) {
//   Matrix<T> augmented = A.augmentedMatrix(b);

//   // Print reduced row echelon form of the augmented matrix for debugging
//   std::cout << "Reduced row echelon form of augmented matrix:" << std::endl;
//   Matrix<T> rrefMatrix = augmented.rref();
//   print_matrix(rrefMatrix);

//   // Check if all pivots are non-zero to determine exact or infinite solution
//   bool allNonZeroPivots = true;
//   for (size_t i = 0; i < rrefMatrix.nRows(); ++i) {
//     if (rrefMatrix(i, i) == 0) {
//       allNonZeroPivots = false;
//       break;
//     }
//   }

//   if (allNonZeroPivots) {
//     // If all pivots are non-zero, calculate exact solution (Xp) using back substitution
//     std::vector<T> xp = backSubstitution(rrefMatrix);
//     return {SolutionType::EXACT_SOLUTION_XP, xp, ""};
//   } else {
//     // If one or more pivots are zero, indicating infinite solutions, Xp should have pivot variables
//     size_t rank = rrefMatrix.Rank();
//     size_t numCols = rrefMatrix.nCols() - 1;
//     std::vector<T> xp(numCols, 0);

//     // Set the pivot variables of Xp to their values from the RREF matrix
//     for (size_t i = 0; i < rank; ++i) {
//       xp[i] = rrefMatrix(i, numCols);
//     }

//     // Calculate Xs for each combination of free variables
//     std::vector<std::string> xsExpressions;
//     size_t numFreeVars = numCols - rank;
//     size_t numXs =
//         numFreeVars == 0 ? 1 : (1 << numFreeVars);  // Number of Xs solutions
//     if (numFreeVars == 1) {
//       numXs =
//           1;  // If there's only one free variable, there should be only one Xs solution
//     }
//     for (size_t xsIdx = 0; xsIdx < numXs; ++xsIdx) {
//       std::vector<T> xs(numCols, 0);
//       size_t freeVarIndex = 0;
//       for (size_t colIdx = rank; colIdx < numCols; ++colIdx) {
//         if ((xsIdx >> freeVarIndex) & 1) {
//            // Set the free variable to 1 if corresponding bit is set
//           xs[colIdx] = 1;
//         }
//         ++freeVarIndex;
//       }

//       // Find the sign of the free variable for each pivot variable
//       std::vector<T> freeVarSigns(rank, 0);  // Initialize all signs to 0
//       for (size_t i = 0; i < rank; ++i) {
//         for (size_t j = 0; j < numCols - 1; ++j) {
//           if (rrefMatrix(i, i) != 0) {
//             freeVarSigns[i] = rrefMatrix(i, j) > 0 ? -1 : 1;
//             break;
//           }
//         }
//       }

//       // Adjust pivot variables according to Xp values
//       for (size_t i = 0; i < rank; ++i) {
//         T pivotValue = rrefMatrix(
//             i,
//             numCols);  // Accessing the last column (b matrix) of the reduced row echelon form matrix
//         T freeVarSign = freeVarSigns[i];
//         std::cout << "Pivot Value: " << pivotValue
//                   << ", Free Variable Sign: " << freeVarSign << std::endl;
//         xs[i] = pivotValue -(freeVarSign * 1);
//       }
//       // If only one free variable, set it to 1 in Xs
//       if (numFreeVars == 1) {
//         xs[numCols - 1] = 1;
//       } else {
//         // If more than one free variable, set one free variable to 1 and others to 0 in Xs
//         size_t freeVarIndex = rank;  // Start from the first free variable index
//         for (size_t i = 0; i < numCols - rank; ++i) {
//           if ((xsIdx >> i) & 1) {
//             // Set the free variable to 1 if corresponding bit is set
//             xs[freeVarIndex++] = 1;
//           } else {
//             // Set the free variable to 0 if corresponding bit is not set
//             xs[freeVarIndex++] = 0;
//           }
//         }
//       }

//       // Format Xs with limited decimal places
//       std::ostringstream xsStream;
//       xsStream << std::fixed << std::setprecision(2);  // Set precision to two decimal places
//       xsStream << "Xs" << xsIdx + 1 << ": ";
//       for (size_t i = 0; i < numCols; ++i) {
//         xsStream << "X" << i + 1 << " = " << xs[i];
//         if (i < numCols - 1)
//           xsStream << "; ";
//       }
//       xsExpressions.push_back(xsStream.str());
//     }

//     // Prepare linear combination of Xs expressions
//     std::string linearComb;
//     for (size_t i = 0; i < xsExpressions.size(); ++i) {
//       linearComb +=
//           "C" + std::to_string(i + 1) + " * (" + xsExpressions[i] + ")";
//       if (i < xsExpressions.size() - 1)
//         linearComb += " + ";
//     }

//     // Return the solution
//     return {SolutionType::INFINITE_SOLUTIONS_XP_XS, xp, linearComb};
//   }
// }

// Function to check if at least one pivot column is zero

// // Define a template function to solve the equation AX = B
// template <typename T>
// std::vector<T> solve_AX_b(const Matrix<T>& A, const std::vector<T>& b) {
//     std::vector<T> solution;

//     try {
//         // Create the augmented matrix [A|B]
//         Matrix<T> augmented = A.augmentedMatrix(b);

//         // Perform row echelonization
//         Matrix<T> echelonMatrix = augmented.rowEchelon();

//         // Check if the system is consistent (has a solution)
//         for (size_t i = 0; i < echelonMatrix.nRows(); ++i) {
//             bool allZero = true;
//             for (size_t j = 0; j < echelonMatrix.nCols() - 1; ++j) {
//                 if (echelonMatrix(i, j) != 0) {
//                     allZero = false;
//                     break;
//                 }
//             }
//             if (allZero && echelonMatrix(i, echelonMatrix.nCols() - 1) != 0) {
//                 // Inconsistent system
//                 throw std::runtime_error("The system is inconsistent and has no solution.");
//             }
//         }

//         // Initialize the solution vector with zeros
//         solution.resize(echelonMatrix.nCols() - 1, T{0});

//         // Backsubstitution to find the solution
//         // Start from the last row and work upwards
//         for (int i = echelonMatrix.nRows() - 1; i >= 0; --i) {
//             T val = echelonMatrix(i, echelonMatrix.nCols() - 1);
//             for (size_t j = i + 1; j < echelonMatrix.nCols() - 1; ++j) {
//                 val -= echelonMatrix(i, j) * solution[j];
//             }
//             solution[i] = val / echelonMatrix(i, i);
//         }
//     } catch (const std::exception& e) {
//         // Catch any exceptions and print an error message
//         std::cerr << "Error: " << e.what() << std::endl;
//     }

//     return solution;
// }

template <typename T>
Matrix<T> REF(const Matrix<T>& input) {
  Matrix<T> result(input);  // Create a copy of the current matrix
  size_t processingRow = 0;

  for (size_t n = 0; n < result.nCols(); ++n) {
    // Check if any entry of the pivot column is zero
    bool pivotColumnZero = true;
    for (size_t m = processingRow; m < result.nRows(); ++m) {
      if (result(m, n) != 0) {
        pivotColumnZero = false;
        break;
      }
    }

    if (pivotColumnZero) {
      continue;  // Skip if pivot column has no nonzero entries
    }

    // Proceed with echelonization process for nonzero pivot column
    size_t pivotRow = processingRow;
    T pivotElement = result(processingRow, n);

    if (pivotElement == 0) {
      // If pivot element is zero, find a nonzero pivot below
      size_t newPivotRow = processingRow + 1;
      while (newPivotRow < result.nRows() && result(newPivotRow, n) == 0) {
        ++newPivotRow;
      }
      if (newPivotRow == result.nRows()) {
        continue;  // No nonzero pivot found, move to the next column
      }
      result.swapRows(processingRow, newPivotRow);
      pivotElement = result(processingRow, n);
    }

    // Perform row operations to eliminate entries below the pivot
    for (size_t q = processingRow + 1; q < result.nRows(); ++q) {
      T localPivotElement = result(q, n);
      for (size_t i = n; i < result.nCols(); ++i) {
        result(q, i) -=
            result(processingRow, i) * localPivotElement / pivotElement;
      }
    }

    processingRow++;
  }

  return result;
}

template <typename T>
Matrix<T> RREF(const Matrix<T>& input) {
  Matrix<T> result(input);  // Create a copy of the current matrix
  size_t processingRow = 0;

  for (size_t col = 0; col < result.nCols(); ++col) {
    // Check if any entry of the pivot column is zero
    bool pivotColumnZero = true;
    for (size_t row = processingRow; row < result.nRows(); ++row) {
      if (result(row, col) != 0) {
        pivotColumnZero = false;
        break;
      }
    }

    if (pivotColumnZero) {
      continue;  // Skip if pivot column has no nonzero entries
    }

    // Proceed with echelonization process for nonzero pivot column
    size_t pivotRow = processingRow;
    T pivotElement = result(processingRow, col);

    if (pivotElement == 0) {
      // If pivot element is zero, find a nonzero pivot below
      size_t newPivotRow = processingRow + 1;
      while (newPivotRow < result.nRows() && result(newPivotRow, col) == 0) {
        ++newPivotRow;
      }
      if (newPivotRow == result.nRows()) {
        continue;  // No nonzero pivot found, move to the next column
      }
      result.swapRows(processingRow, newPivotRow);
      pivotElement = result(processingRow, col);
    }

    // Normalize the pivot row to make the pivot element 1
    for (size_t i = col; i < result.nCols(); ++i) {
      if (result(processingRow, i) != 0 && pivotElement != 0) {
        result(processingRow, i) /= pivotElement;
      }
    }

    // Perform row operations to eliminate entries above and below the pivot
    for (size_t q = 0; q < result.nRows(); ++q) {
      if (q != processingRow && result(q, col) != 0) {
        T localPivotElement = result(q, col);
        for (size_t i = col; i < result.nCols(); ++i) {
          result(q, i) -= result(processingRow, i) * localPivotElement;
        }
      }
    }

    processingRow++;
  }

  return result;
}

#undef EPSILON
#undef EQUAL
#undef INDEX
#undef for_size
#undef for_ij
#undef CHECK_SIZE
#undef CHECK_ROW
#undef CHECK_COLUMN
#undef VALIDATE_INDEX

#endif
