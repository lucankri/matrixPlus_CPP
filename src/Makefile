GCC = g++ -Wall -Werror -Wextra

all: clean s21_matrix_oop.a test

.PHONY: all test s21_matrix_oop.a gcov_report

clean:
	rm -rf s21_matrix_oop.o s21_matrix_oop.a test
	rm -rf *.gcda *.gcno report *.info
	rm -rf .clang-format

s21_matrix_oop.a:
	$(GCC) -c s21_matrix_oop.cpp -o s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

test: s21_matrix_oop.a
	$(GCC) tests.cpp -lgtest s21_matrix_oop.a -o test
	./test

gcov_report: s21_matrix_oop.a
	$(GCC) --coverage -c s21_matrix_oop.cpp -o s21_matrix_oop.o
	$(GCC) --coverage tests.cpp -lgtest s21_matrix_oop.o -o test
	./test
	lcov -t "tests" -o tests.info -c -d s21_matrix_oop.gcda
	genhtml -o report tests.info

check:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -n *.cpp
	clang-format -n *.h
	cppcheck --language=c++ *.cpp
	cppcheck --language=c++ *.h
