#include "s21_matrix_oop.h"

/*--------------------Конструкторы--------------------*/

S21Matrix::S21Matrix() noexcept {  // Базовый конструтор
  this->rows_ = 0;
  this->cols_ = 0;
  this->matrix_ = nullptr;
}

S21Matrix::S21Matrix(int rows_, int cols_) {  // Параметризированный конструктор
  this->create(rows_, cols_);
}

S21Matrix::S21Matrix(const S21Matrix& other) {  // Конструктор копирования
  this->create(other.rows_, other.cols_);
  for (int i = 0; i < other.rows_; i++) {
    std::memmove(this->matrix_[i], other.matrix_[i],
                 sizeof(double) * other.cols_);
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept {  // Конструктор перемещения
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {  // Деструктор
  this->clear();
}

/*--------------------Методы--------------------*/

/**
 * @brief Проверяет матрицы на равенство между собой
 *
 * @param other
 * @return true
 * @return false
 */
bool S21Matrix::EqMatrix(const S21Matrix& other) const noexcept {
  bool tmp = this->rows_ == other.rows_ && this->cols_ == other.cols_;
  for (int i = 0; i < this->rows_ && tmp == true; i++) {
    for (int j = 0; j < this->cols_ && tmp == true; j++) {
      if (std::fabs(this->matrix_[i][j] - other.matrix_[i][j]) >
          S21Matrix::EPS) {
        tmp = false;
      }
    }
  }
  return tmp;
}

/**
 * @brief Прибавляет вторую матрицу к текущей
 *
 * @param other
 */
void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_)
    throw Exception("Different dimensions of matrices!");
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

/**
 * @brief Вычитает из текущей матрицы другую
 *
 * @param other
 */
void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_)
    throw Exception("Different dimensions of matrices!");
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

/**
 * @brief Умножает текущую матрицу на число
 *
 * @param num
 */
void S21Matrix::MulNumber(const double num) noexcept {
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->cols_ != other.rows_)
    throw Exception(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix!");
  S21Matrix tmp(this->rows_, other.cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      tmp.matrix_[i][j] = 0;
      for (int k = 0; k < this->cols_; k++) {
        tmp.matrix_[i][j] += this->matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = std::move(tmp);
}

/**
 * @brief Создает новую транспонированную матрицу из текущей и возвращает ее
 *
 * @return S21Matrix
 */
S21Matrix S21Matrix::Transpose() const {
  S21Matrix result(this->cols_, this->rows_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      result.matrix_[j][i] = this->matrix_[i][j];
    }
  }
  return result;
}

/**
 * @brief Вычисляет матрицу алгебраических дополнений текущей матрицы и
 * возвращает ее
 *
 * @return S21Matrix
 */
S21Matrix S21Matrix::CalcComplements() const {
  if (this->rows_ != this->cols_) throw Exception("Matrix is not square!");
  S21Matrix result(this->rows_, this->cols_);
  if (this->rows_ != 1) {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        S21Matrix tmp = this->minor(i, j);
        result.matrix_[i][j] = tmp.Determinant() * std::pow(-1, i + j);
      }
    }
  } else {
    result.matrix_[0][0] = 1;
  }
  return result;
}

/**
 * @brief Вычисляет и возвращает определитель текущей матрицы
 *
 * @return double
 */
double S21Matrix::Determinant() const {
  if (this->rows_ != this->cols_) throw Exception("Matrix is not square!");
  /*Алгоритм треугольной матрицы с перестановками...*/

  S21Matrix M = *this;
  double result = 1;
  int N = M.rows_;
  int exchanges = 0;
  double denominator = 1;
  for (int i = 0; i < N - 1 && result == 1; i++) {
    int maxN = i;
    double maxValue = std::fabs(M.matrix_[i][i]);
    for (int j = i + 1; j < N; j++) {
      double value = std::fabs(M.matrix_[j][i]);
      if (value > maxValue) {
        maxN = j;
        maxValue = value;
      }
    }
    if (maxN > i) {
      double* tmp = M.matrix_[i];
      M.matrix_[i] = M.matrix_[maxN];
      M.matrix_[maxN] = tmp;
      exchanges++;
    } else {
      if (maxValue == 0) {
        result = maxValue;
      }
    }
    if (result == 1) {
      double valueAbove = M.matrix_[i][i];
      for (int j = i + 1; j < N; j++) {
        double valueLeft = M.matrix_[j][i];
        M.matrix_[j][i] = 0;
        for (int k = i + 1; k < N; k++)
          M.matrix_[j][k] =
              (M.matrix_[j][k] * valueAbove - M.matrix_[i][k] * valueLeft) /
              denominator;
      }
      denominator = valueAbove;
    }
  }
  if (result == 1 && M.matrix_[N - 1][N - 1] != 0) {
    exchanges % 2 ? result = -M.matrix_[N - 1][N - 1]
                  : result = M.matrix_[N - 1][N - 1];
  } else {
    result = 0;
  }
  return result;
}

/**
 * @brief Вычисляет и возвращает обратную матрицу
 *
 * @return S21Matrix
 */

S21Matrix S21Matrix::InverseMatrix() const {
  double determinant = this->Determinant();
  if (determinant == 0) throw Exception("Inverse matrix does not exist!");

  S21Matrix result = this->CalcComplements();
  result = result.Transpose();
  for (int i = 0; i < result.rows_; i++) {
    for (int j = 0; j < result.cols_; j++) {
      result.matrix_[i][j] /= determinant;
    }
  }
  return result;
}

/*--------------------Операторы--------------------*/

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) return *this;

  this->clear();
  this->create(other.rows_, other.cols_);
  for (int i = 0; i < other.rows_; i++) {
    std::memmove(this->matrix_[i], other.matrix_[i],
                 sizeof(double) * other.cols_);
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
  if (this == &other) return *this;

  this->clear();
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;

  return *this;
}

double& S21Matrix::operator()(const int& i, const int& j) const {
  if (!(i < this->rows_ && j < this->cols_) || (i < 0 || j < 0))
    throw Exception("Out of bounds array!");
  return this->matrix_[i][j];
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix result = *this;
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix result = *this;
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix result = *this;
  result.MulMatrix(other);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return this->EqMatrix(other);
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double value) {
  this->MulNumber(value);
  return *this;
}

/*--------------------Оперторы вне класса--------------------*/

S21Matrix operator*(S21Matrix& Matrix, double value) {
  S21Matrix result = Matrix;
  result.MulNumber(value);
  return result;
}

S21Matrix operator*(double value, S21Matrix& Matrix) {
  S21Matrix result = Matrix * value;
  return result;
}

/*--------------------"Геттеры" и "Сеттеры"--------------------*/

/**
 * @brief Метод выводит количество строк
 *
 * @return int
 */
int S21Matrix::GetRows() const noexcept { return this->rows_; }

/**
 * @brief Метод выводит количество колонок
 *
 * @return int
 */
int S21Matrix::GetCols() const noexcept { return this->cols_; }

/**
 * @brief Метод увеличивает количество строк в матрице.
 * Матрица дополняется нулевыми элементами
 * При увеличении размера - матрица дополняется нулевыми элементами, при
 * уменьшении - лишнее просто отбрасывается
 *
 * @param rows_
 */
void S21Matrix::SetRows(int rows_) {
  S21Matrix result(rows_, this->cols_);
  for (int i = 0; i < rows_; i++) {
    if (i < this->rows_) {
      std::memmove(result.matrix_[i], this->matrix_[i],
                   sizeof(double) * this->cols_);
    } else {
      std::memset(result.matrix_[i], 0, sizeof(double) * this->cols_);
    }
  }
  *this = std::move(result);
}

/**
 * @brief Метод увеличивыет или уменьшает количество колонок в матрице.
 * При увеличении размера - матрица дополняется нулевыми элементами, при
 * уменьшении - лишнее просто отбрасывается
 *
 * @param cols_
 */
void S21Matrix::SetCols(int cols_) {
  S21Matrix result(this->rows_, cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (j < this->cols_) {
        result.matrix_[i][j] = this->matrix_[i][j];
      } else {
        result.matrix_[i][j] = 0;
      }
    }
  }
  *this = std::move(result);
}

/*--------------------Методы private--------------------*/

/**
 * @brief private Метод найти минор матрицы по rows и cols
 *
 * @param rows_
 * @param cols_
 * @return S21Matrix
 */
S21Matrix S21Matrix::minor(int row_minor, int col_minor) const {
  S21Matrix tmp(this->rows_ - 1, this->cols_ - 1);
  for (int i = 0, row = 0; i < this->rows_; i++) {
    if (i != row_minor) {
      for (int j = 0, col = 0; j < this->cols_; j++) {
        if (j != col_minor) {
          tmp.matrix_[row][col++] = this->matrix_[i][j];
        }
      }
      row++;
    }
  }
  return tmp;
}

/**
 * @brief private Метод создание матрицы;
 *
 * @param rows
 * @param cols
 */
void S21Matrix::create(int rows, int cols) {
  this->rows_ = rows;
  this->cols_ = cols;
  if (rows == 0 && cols == 0) {
    this->rows_ = 0;
    this->cols_ = 0;
    matrix_ = nullptr;
  } else if (rows == 0 || cols == 0) {
    throw Exception("One of the parameters 0. Cannot allocate memory!");
  } else {
    try {
      matrix_ = new double*[rows_];
      for (int i = 0; i < rows; i++) {
        try {
          matrix_[i] = new double[cols];
        } catch (...) {
          for (int j = 0; j < i; j++) {
            delete[] matrix_[i];
          }
          delete[] matrix_;
          throw Exception("Memory is not allocated!");
        }
      }
    } catch (...) {
      throw Exception("Memory is not allocated!");
    }
  }
}

/**
 * @brief private Метод очистки матрицы;
 */
void S21Matrix::clear() {
  if (this->matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) delete[] this->matrix_[i];
    delete[] this->matrix_;
  }
}
