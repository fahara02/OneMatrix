// //Blocked_Bidiagonal_Reduction_Howell
// Matrix<T> BidiagRed_Howell(size_t block_size) {

//   Matrix<T> A = *this;
//   size_t m = A.nRows();
//   size_t n = A.nCols();
//   size_t current = 0;

//   Matrix<T> U, V, Y, Z;

//   while (current < std::min(m, n)) {
//     size_t b = std::min(block_size, std::min(m, n) - current);

//     Matrix<T> A_block = A.slice(current, m, current, n);
//     Matrix<T> U_block, V_block, Y_block, Z_block;

//     BiRed_Lazy_Rearranged_Unblocked(A_block, U_block, V_block, Y_block,
//                                     Z_block, b);

//     // Update trailing submatrix
//     if (current + b < std::min(m, n)) {
//       Matrix<T> A22 = A.slice(current + b, m, current + b, n);
//       Matrix<T> U2 = U_block.slice(b, U_block.nRows(), 0, U_block.nCols());
//       Matrix<T> Y2 = Y_block.slice(b, Y_block.nRows(), 0, Y_block.nCols());
//       Matrix<T> Z2 = Z_block.slice(b, Z_block.nRows(), 0, Z_block.nCols());
//       Matrix<T> V2 = V_block.slice(b, V_block.nRows(), 0, V_block.nCols());

//       A22 = A22 - U2 * Y2.transpose() - Z2 * V2.transpose();
//     }

//     current += b;
//   }
//   return A;
// }

// void BiRed_Lazy_Rearranged_Unblocked(Matrix<T>& A, Matrix<T>& U, Matrix<T>& V,
//                                      Matrix<T>& Y, Matrix<T>& Z, size_t b) {
//   size_t m = A.nRows();
//   size_t n = A.nCols();

//   U = Matrix<T>::zeros(m, b);
//   V = Matrix<T>::zeros(n, b);
//   Y = Matrix<T>::zeros(m, b);
//   Z = Matrix<T>::zeros(n, b);

//   for (size_t k = 0; k < b; ++k) {
//     // Extract current column and row
//     T alpha = A(k, k);
//     Matrix<T> a_col = A.slice(k + 1, m, k, k + 1);
//     Matrix<T> a_row = A.slice(k, k + 1, k + 1, n);

//     // Update with previous transformations
//     if (k > 0) {
//       Matrix<T> U_prev = U.slice(k, m, 0, k);
//       Matrix<T> Y_prev = Y.slice(k, m, 0, k);
//       Matrix<T> Z_prev = Z.slice(k, m, 0, k);
//       Matrix<T> V_prev = V.slice(k + 1, n, 0, k);

//       a_col = a_col - U_prev * Y_prev.transpose().getColumn(k) -
//               Z_prev * V_prev.transpose().getColumn(k);
//       a_row = a_row -
//               U.slice(k, k + 1, 0, k).transpose() * Y.slice(k + 1, n, 0, k) -
//               Z.slice(k, k + 1, 0, k).transpose() * V.slice(k + 1, n, 0, k);
//     }

//     // Compute Householder vector for column
//     Matrix<T> col_vec = Matrix<T>::vstack(Matrix<T>(1, 1, alpha), a_col);
//     auto [u_col, tau_col] = householderColumn(col_vec);
//     U.setColumn(u_col.slice(0, m, 0, 1), k);
//     alpha = u_col(0, 0);
//     A(k, k) = alpha;
//     A.slice(k + 1, m, k, k + 1) = u_col.slice(1, u_col.nRows(), 0, 1);

//     // Compute Householder vector for row
//     Matrix<T> row_vec = Matrix<T>::hstack(Matrix<T>(1, 1, alpha), a_row);
//     auto [v_row, tau_row] = householderRow(row_vec);
//     V.setColumn(v_row.slice(0, n, 0, 1), k);
//     A.slice(k, k + 1, k + 1, n) = v_row.slice(1, v_row.nCols(), 0, 1);

//     // Compute Y and Z updates
//     Y.setColumn((A.transpose() * U.getColumn(k)) / tau_col, k);
//     Z.setColumn((A * V.getColumn(k) - alpha * A.getColumn(0)) / tau_row, k);
//   }
// }

// void Lazy_Unblocked_Bidiag_Panel(Matrix<T>& A_panel, Matrix<T>& A_row,
//                                  size_t nb, Matrix<T>& U_panel,
//                                  Matrix<T>& V_panel, Matrix<T>& T_L,
//                                  Matrix<T>& T_R) {
//   size_t m = A_panel.nRows();
//   size_t n_row = A_row.nCols();

//   U_panel = Matrix<T>::zeros(m, nb);
//   V_panel = Matrix<T>::zeros(n_row, nb);
//   T_L = Matrix<T>::zeros(nb, nb);
//   T_R = Matrix<T>::zeros(nb, nb);

//   for (size_t j = 0; j < nb; ++j) {
//     // Left transformation
//     Matrix<T> x = A_panel.slice(j, m, j, j + 1);
//     auto [u_j, tau_L, x_new] = Housev(x);

//     // Update A_panel column
//     for (size_t i = j; i < m; ++i)
//       A_panel(i, j) = x_new(i - j, 0);

//     // Store Householder vector
//     for (size_t i = j; i < m; ++i)
//       U_panel(i, j) = u_j(i - j, 0);

//     T_L(j, j) = tau_L;

//     // Apply left to trailing columns
//     if (j < nb - 1) {
//       Matrix<T> trailing = A_panel.slice(j, m, j + 1, nb);
//       Matrix<T> update = (u_j.transpose() * trailing) * (1 / tau_L);
//       trailing = trailing - u_j * update;
//       // Update A_panel
//       for (size_t r = 0; r < trailing.nRows(); ++r)
//         for (size_t c = 0; c < trailing.nCols(); ++c)
//           A_panel(j + r, j + 1 + c) = trailing(r, c);
//     }

//     // Right transformation
//     if (j < n_row - 1) {
//       Matrix<T> row = A_row.slice(j, j + 1, j + 1, n_row);
//       auto [v_j, tau_R, a_new] = Housev(row.transpose());

//       // Update A_row
//       for (size_t c = 0; c < a_new.nRows(); ++c)
//         A_row(j, j + 1 + c) = a_new(c, 0);

//       // Store Householder vector
//       for (size_t i = j + 1; i < n_row; ++i)
//         V_panel(i, j) = v_j(i - (j + 1), 0);

//       T_R(j, j) = tau_R;

//       // Apply right to trailing rows
//       if (j < nb - 1) {
//         Matrix<T> trailing = A_panel.slice(j + 1, m, j + 1, nb);
//         Matrix<T> update = (trailing * v_j) * (1 / tau_R);
//         trailing = trailing - update * v_j.transpose();
//         // Update A_panel
//         for (size_t r = 0; r < trailing.nRows(); ++r)
//           for (size_t c = 0; c < trailing.nCols(); ++c)
//             A_panel(j + 1 + r, j + 1 + c) = trailing(r, c);
//       }
//     }
//   }

//   // Compute T_L and T_R
//   for (size_t i = 0; i < nb; ++i) {
//     for (size_t j = i + 1; j < nb; ++j) {
//       Matrix<T> ui = U_panel.getColumn(i);
//       Matrix<T> uj = U_panel.getColumn(j);
//       T_L(i, j) = (ui.transpose() * uj)(0, 0) / 2;

//       Matrix<T> vi = V_panel.getColumn(i);
//       Matrix<T> vj = V_panel.getColumn(j);
//       T_R(i, j) = (vi.transpose() * vj)(0, 0) / 2;
//     }
//   }
// }