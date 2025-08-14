// Matrix header

#include <vector>

template<typename T>
class Matrix {
public:
	size_t rows, cols;
	std::vector<T> data;

	Matrix(size_t rows, size_t cols, const T& initVal = T())
        : rows(rows), cols(cols), data(rows* cols, initVal) {}

    T& operator()(size_t r, size_t c) {
        if (r >= rows || c >= cols)
            throw std::out_of_range("Matrix index out of range");
        return data[r * cols + c];
    }

    const T& operator()(size_t r, size_t c) const {
        if (r >= rows || c >= cols)
            throw std::out_of_range("Matrix index out of range");
        return data[r * cols + c];
    }


    Matrix<T> operator+(const Matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for addition");

        Matrix<T> result(rows, cols);
        for (size_t i = 0; i < data.size(); ++i)
            result.data[i] = data[i] + other.data[i];
        return result;
    }

    Matrix<T> operator*(const Matrix<T>& other) const {
        if (cols != other.rows)
            throw std::invalid_argument("Matrix dimensions incompatible for multiplication");

        Matrix<T> result(rows, other.cols, T());
        for (size_t i = 0; i < rows; ++i)
            for (size_t j = 0; j < other.cols; ++j)
                for (size_t k = 0; k < cols; ++k)
                    result(i, j) += (*this)(i, k) * other(k, j);
        return result;
    }

    Matrix<T> operator*(const T& scalar) const {
        Matrix<T> result(rows, cols);
        for (size_t i = 0; i < data.size(); ++i)
            result.data[i] = data[i] * scalar;
        return result;
    }

    Matrix<T> tensorProduct(const Matrix<T>& other) const {
        Matrix<T> result(rows * other.rows, cols * other.cols);

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                T value = (*this)(i, j);
                for (size_t ii = 0; ii < other.rows; ++ii) {
                    for (size_t jj = 0; jj < other.cols; ++jj) {
                        result(i * other.rows + ii, j * other.cols + jj) =
                            value * other(ii, jj);
                    }
                }
            }
        }

        return result;
    }


    void print() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j)
                std::cout << (*this)(i, j) << " ";
            std::cout << "\n";
        }
    }
};