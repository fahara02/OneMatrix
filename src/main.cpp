#include <iomanip>
#include <iostream>
#include <ratio>
#include <sstream>

#include "matrix.h"

int main() {
  // Matrix<int> mat(3, 4, {1, 2, 1, 1, 2, 4, 2, 2, 3, 6, 3, 4});
  // Matrix<double> mat2(3, 3, {3, 1, 0, -2, -4, 3, 5, 4, -2});
  //   Matrix<double> big_mat(
  //       4, 5, {1, 2, 2, 3, 1, 2, 4, 4, 6, 2, 3, 6, 6, 9, 6, 1, 2, 4, 5, 3});

  // Matrix<double> big_mat2(5, 5, {1, 2, 1, 3, 3, 2, 4, 0, 4, 4, 1, 2, 3,
  //                                5, 5, 2, 4, 0, 4, 7, 1, 2, 2, 3, 1});

  // // Matrix<double> big_mat3(3, 5,
  // //                         {2, -4, -8, 6, 3, 0, 1, 3, 2, 3, 3, -2, 0, 0, 8});
  // // std::cout<<"printing the matrix"<<std::endl;
  // // print_matrix(big_mat);
  // // Matrix<int> mat2 = {// Create a 3x2 matrix with initializer list
  // //                     {1, 2},
  // //                     {3, 4},
  // //                     {5, 6}};

  // // // Create a vector containing the data for the matrix
  // // std::vector<int> data = {11, 21, 31, 41, 51, 61};
  // // Matrix<int> mat3(3, 2, data);
  // // Matrix<int> mat4(3, 2, {12, 22, 32, 42, 52, 72});

  // // int data2[6] = {1, 21, 3, 41, 5, 71};  // Array containing the data

  // // Matrix<int> mat5(2, 3, data2);  // Create a 3x2 matrix using the constructor
  // // Matrix<int> bin(3, 1, {1, 0, 1});
  // // Matrix<int> row(1, 3, {1, 0, 1});
  // // Matrix<int> col(3, 1, {1, 0, 1});
  // // int a = 2;
  // // Matrix<int> newRow(1, 3);
  // // newRow(0, 0) = 3;
  // // newRow(0, 1) = 4;
  // // newRow(0, 2) = 5;
  // // mat.addRow(newRow, 2);
  // // mat += 2;
  // // mat *= 12;
  // // std::cout << mat.Size();

  // // print_matrix(mat.getRows(bin));
  // //  print_matrix(mat.getColumns(bin));
  // // mat.setRow(0, row);

  // // print_matrix(mat2);
  // // mat.setColumn(1, col);

  // // mat.removeColumn(0);

  // // mat.swapColumns(1, 2);
  // // std::cout << "index of :" << get<0>(mat.min()) << "," << get<1>(mat.min())
  // //           << " is minimum :=" << get<2>(mat.min());
  // // Matrix<int> rows = mat.getRow(0);
  // // std::cout << "printing original matrix" << std::endl;
  // // print_matrix(big_mat);

  // // std::cout << "rows no=" << rows.nRows() << std::endl;
  // // std::cout << "cols no=" << rows.nCols();

  // // Matrix<int> cols = mat.getColumn(0);
  // // print_matrix(cols);
  // // std::cout << "rows no=" << cols.nRows() << std::endl;
  // // std::cout << "cols no=" << cols.nCols();
  // // std::cout << "printing swap matrix" << std::endl;

  // // print_matrix(E01);
  // // std::cout << std::endl;

  // // Matrix<int> mate = mat.echelon(true);
  // // std::cout << "printing echelon" << std::endl;
  // // print_matrix(mate);

  // // Matrix<int> mat3 = mat.echelon(true);
  // // std::cout << "printing echelon matrix" << std::endl;
  // // print_matrix(mat3);
  // // std::vector<Matrix<double>> rowOps;
  // // Matrix<double> mat3 = REF(mat2);
  // // std::cout << " row echelon form is" << std::endl;
  // // print_matrix(mat3);

  // // std::vector<Matrix<double>> rowOperations;

  // // Matrix<double> mat4 = mat2.echelon_form();
  // //  std::cout << "finally echelon matrix is" << std::endl;
  // // print_matrix(mat4);

  // // std::cout << "rank is " << mat4.Rank() << std::endl;

  // // print_matrix(mat2.identity());
  // // print_matrix(mat2.zeros());
  // // std::cout << " row echelon from original function is " << std::endl;
  // //
  // // std::cout << " row echelon from new broken down function  is" << std::endl;
  // // print_matrix(big_mat.row_echelon_form());

  // // std::cout << " reduced row echelon from original function is " << std::endl;
  // // print_matrix(big_mat.echelon_big(true));
  // // std::cout << " reduced row echelon from new broken down function  is"
  // //           << std::endl;
  // // print_matrix(big_mat.reduced_row_echelon_form());

  // // //std::cout << "reduced row echelon matrix is" << std::endl;
  // // // print_matrix(big_mat.reduced_echelon_form());
  // // std::cout << "rank is:" << big_mat.Rank() << std::endl;
  // //  print_matrix(big_mat.rowEchelon());
  // //  print_matrix(big_mat.rref());
  // //  std::cout<<"the determinant of the matrix is :"<< mat2.det()<<std::endl;

  // Matrix<double> v1({1, 0, 0}, VectorType::RowVector);
  //  Matrix<double> v2({3, 1, 2}, VectorType::RowVector);

  //  Matrix<double> v3({1, 0, -2}, VectorType::RowVector);
  // Matrix<double> v4({2, 0, 0}, VectorType::ColumnVector);

  //  Matrix<double> C({v1, v2, v3}, MatrixOrientation::Row);

  //   std::cout << "===========main matrix===============" << std::endl;
  //   print_matrix(C);

  //   std::cout << "==========================" << std::endl;

  //   std::cout << "printing the  matrix " << std::endl;
  //   print_matrix(C);
  //   std::cout << "printing the echelon matrix " << std::endl;
  //   print_matrix(C.rowEchelon());

  // // Call the function to get pivot columns
  // std::vector<std::pair<size_t, Matrix<int>>> pivotColumns = C.getPivotColumns();

  // // Output the results
  // for (const auto& pivotColumn : pivotColumns) {
  // 	std::cout << "Pivot Column Index: " << pivotColumn.first << std::endl;
  // 	print_matrix(pivotColumn.second);
  // }
  // std::cout << "the rank of the matrix is :" << C.Rank() << std::endl;
  // std::cout << "the determinant of the matrix is :" << C.det() << std::endl;

  // std::cout << "printing a new matrix" << std::endl;
  // print_matrix(mat);
  // Matrix<int> vColumn({1, 0, -1}, VectorType::ColumnVector);
  // mat.addColumn(vColumn, 2);
  // mat.addColumn({2, 0, 0, 3});
  // std::cout << "printing after row addition" << std::endl;
  // print_matrix(mat);
  // std::cout << "rows=" << mat.nRows() << "columns=" << mat.nCols() << std::endl;
  // std::cout << "printing swapped rows" << std::endl;
  // print_matrix(mat.swapRows(1, 2));

  // std::cout << "printing swapped column" << std::endl;
  // print_matrix(mat.swapColumns(1, 2));

  // try {
  // 	auto val1 = mat.get_element(1, 2);
  // 	std::cout << val1 << std::endl;
  // 	// Use 'value' here
  // }
  // catch (const std::exception& e) {
  // 	std::cerr << "Error: " << e.what() << std::endl;
  // 	// Handle the error gracefully, e.g., return a default value, throw another
  // 	// exception, or terminate the program
  // }

  // Fill the matrix with values

  // Fill the matrix with values

  // Matrix<int> mat(3, 3);
  // try {
  //   size_t value = -1;
  //   for (auto it = mat.begin(); it != mat.end(); ++it) {
  //     *it = ++value;
  //     //std::cout << "Index: " << *it << std::endl;
  //   }

  // std::cout << "=========== matrix ===============" << std::endl;
  // print_matrix(big_mat);

  // std::cout << "=============original rref=============" << std::endl;
  // print_matrix(big_mat.rref());

  // } catch (const std::exception& e) {
  //   std::cerr << "Error: " << e.what() << std::endl;
  // }

  // for (auto it = --mat.end(); it != mat.begin();) {
  //   --it;  // Decrement before using the iterator
  //   std::cout << "row" << it.getrow() << " col" << it.getcol() << "=" << *it
  //             << std::endl;
  // }

  // std::cout << "printing all elements again" << std::endl;
  // for (auto& element : mat) {
  //   std::cout << element << " ";
  // }
  // // Advancing the iterator by a certain number of steps
  // auto it = mat.begin();
  // it.advance(2);  // Advance by 5 elements
  // std::cout << "Value after advancing: " << *it << std::endl;

  // // Calculating distance between iterators
  // auto it1 = mat.begin();
  // auto it2 = mat.begin() + 2;
  // std::cout << "Distance between iterators: " << it1.distance_to(it2) << std::endl;

  // Suppose you want to iterate over a submatrix of the original matrix
  // Create a submatrix iterator
  // size_t startRow = 1;
  // size_t startCol = 1;
  // size_t numRows = 2;
  // size_t numCols = 2;
  // auto submat_begin = mat.begin() + startRow * mat.nCols() + startCol;
  // auto submat_end = submat_begin + numRows * mat.nCols() + numCols;
  // std::cout << "all subs" << std::endl;
  // // Iterate over the submatrix
  // for (auto it = submat_begin; it != submat_end; ++it) {

  // 	std::cout << *it << " ,";
  // }

  //  auto [L,D, U] = C.to_LDU();

  // Print U

  //     // Print L
  //     std::cout << "L matrix:" << std::endl;
  //  print_matrix(L);
  // //     std::cout << "D matrix:" << std::endl;
  // //  print_matrix(D);

  //     std::cout << "output U matrix:" << std::endl;
  //  print_matrix(U);

  //     std::cout << "D" << std::endl;
  //  print_matrix(D);
  //        // Print L
  //     std::cout << " from L*D*U " << std::endl;
  //  print_matrix(L*D*U);

  //  std::cout << " determinant is " << C.det() << std::endl;

  // Matrix<double> Inv= C.augmentedMatrix(C.to_identity());
  //  std::cout << " original matrix "<< std::endl;
  //   print_matrix(Inv);
  //     std::cout << " ----rref matrix-----" << std::endl;
  //   print_matrix(Inv.rref());
  //   Matrix<double> bin({0, 0, 0,1,1,1}, VectorType::ColumnVector);
  //   Matrix<double> Inverse=Inv.rref().getColumns(bin);
  //   std::cout << " ----inverse matrix-----" << std::endl;
  //    print_matrix(Inverse);

  // Matrix<float> floatMatrix(3, 3, {1.5233333333333333333f, 2.5f, 3.5f, 4.5f, 5.5f, 6.5f, 7.5f, 8.5f, 9.5f});

  // Matrix<int> v1({1, 0, 0}, VectorType::RowVector);
  // Matrix<int> v2({3, 1, 2}, VectorType::RowVector);
  // Matrix<int> v3({1, 0, -2}, VectorType::RowVector);

  // Matrix<int> C({v1, v2, v3}, Orientation::Row);
  // print_matrix(C);
  // print_matrix(floatMatrix);

  //  std::cout << " ----addition-----" << std::endl;
  // // floatMatrix += C;
  // // print_matrix(floatMatrix);

  // C+=floatMatrix;

  // print_matrix(C);
  // Test matrix and vector
  // Test matrix A

  // Print the matrix A

  // Solve the matrix equation

  // try{  Solution<double> solution = solve_AX_b(A, b);
  // // Print the solution
  // std::cout << "\nSolution:\n";
  // if (solution.size() == 1 && solution[0] == 0) {
  //   std::cout << "No solution\n";
  // } else {
  //   for (size_t i = 0; i < solution.size(); ++i) {
  //     std::cout << "x" << i + 1 << " = " << solution[i] << "\n";
  //   }
  // }  }catch (const std::exception& e) {
  //       // Catch any exceptions and print an error message
  //       std::cerr << "Error: " << e.what() << std::endl;
  //   }

  //  // Solve the equation AX = b
  //   Solution<double> solution = solve_AX_b(A, b);

  //   // Print the solution
  //   switch (solution.type) {
  //       case SolutionType::EXACT_SOLUTION:
  //           std::cout << "Exact solution:\n";
  //           for (size_t i = 0; i < solution.values.size(); ++i) {
  //               std::cout << "x" << i + 1 << " = " << solution.values[i] << "\n";
  //           }
  //           break;
  //       case SolutionType::INFINITE_SOLUTIONS:
  //           std::cout << "Infinite solutions." << std::endl;
  //           break;
  //       case SolutionType::NO_SOLUTION:
  //           std::cout << "No solution found." << std::endl;
  //           break;
  //   }
  float Kr = 0.299f, Kg = 0.587f, Kb = 0.114f;
  float uScale = 0.5f, vScale = 0.5f;

  // Initialize using nested initializer lists
  Matrix<float> rgbToYUV = {
      {Kr, Kg, Kb},
      {-Kr / (1.0f - Kb) * uScale, -Kg / (1.0f - Kb) * uScale, uScale},
      {vScale, -Kg / (1.0f - Kr) * vScale, -Kb / (1.0f - Kr) * vScale}};

  // Print the matrix to verify
  print_matrix(rgbToYUV);

  // Define matrix A and vector b
  Matrix<double> A(3, 4, {1, 3, 0, 2, 0, 0, 1, 4, 1, 3, 1, 6});
  // // Matrix<double> A(3, 3, {1, 1, 2,
  //                       2, 4, -3,
  //                       3, 6, -5});
  std::cout << "Original matrix:\n";
  print_matrix(A);

  std::vector<double> b = {1, 6, 7};
  std::cout << "Matrix A augmented:\n";
  print_matrix(A.augmentedMatrix(b));

  std::cout << "Matrix A augmented rref:\n";
  print_matrix(A.augmentedMatrix(b).rref());

  try {
    Solution<double> solution = solver_AX_b(A, b);

    // Printing the solution
    if (solution.type == SolutionType::EXACT_SOLUTION) {
      std::cout << "Exact solution (Xp):\n";
      print_matrix(solution.XP);
    } else if (solution.type == SolutionType::INFINITE_SOLUTIONS) {
      std::cout << "Infinite solutions. Particular solution (Xp):\n";
      print_matrix(solution.XP);
      print_matrix(solution.XS[0]);
      print_matrix(solution.XS[1]);
      std::cout << "Linear combination of special solutions (Xs):\n"
                << solution.linearComb << std::endl;
    } else {
      std::cout << "No solution exists." << std::endl;
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
