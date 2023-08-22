#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_
#include <cmath>
#include <cstring>
#include <exception>
#include <new>
#include <utility>

class Exception : public std::exception {
 public:
  Exception(const char* msg) { this->error = msg; }
  const char* GetError() { return error; }

 private:
  const char* error;
};

class S21Matrix {
  friend S21Matrix operator*(S21Matrix& Matrix, double value);
  friend S21Matrix operator*(double value, S21Matrix& Matrix);

 public:
  /*-----Конструкторы-----*/

  S21Matrix() noexcept;  // Базовый конструтор
  S21Matrix(int rows_, int cols_);  // Параметризированный конструктор
  S21Matrix(const S21Matrix& other);  // Конструктор копирования
  S21Matrix(S21Matrix&& other) noexcept;  // Конструктор перемещения
  ~S21Matrix();                           // Деструктор

  /*-----Операторы-----*/

  S21Matrix& operator=(
      const S21Matrix& other);  // Оператор присваивания копированием
  S21Matrix& operator=(
      S21Matrix&& other) noexcept;  // Оператор присваивания перемещением
  double& operator()(const int& i,
                     const int& j) const;  // Оператор индексации rows и cols
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const S21Matrix& other) const;
  bool operator==(const S21Matrix& other) const;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double value);

  /*-----Методы-----*/

  bool EqMatrix(const S21Matrix& other) const noexcept;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num) noexcept;
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;

  /*-----"Геттеры" и "Сеттеры"-----*/

  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetRows(int rows_);
  void SetCols(int cols_);

 private:
  int rows_, cols_;
  double** matrix_;

  static constexpr double EPS = 1e-7;

  /*-----private Методы-----*/

  S21Matrix minor(int rows_minor, int cols_minor) const;
  void create(int rows, int cols);
  void clear();
};

/*-----Операторы вне класса-----*/

S21Matrix operator*(S21Matrix& Matrix, double value);
S21Matrix operator*(double value, S21Matrix& Matrix);

#endif  // SRC_S21_MATRIX_OOP_H_
