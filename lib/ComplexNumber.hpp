#include <cmath>
#include <iostream>

class ComplexNumber {
private:
    double m_real;
    double m_imag;

public:
    ComplexNumber(double r = 0.0, double i = 0.0) : m_real(r), m_imag(i) {}

    ComplexNumber operator+(const ComplexNumber& other) const {
        return ComplexNumber(m_real + other.m_real, m_imag + other.m_imag);
    }

    ComplexNumber& operator+=(const ComplexNumber& other) {
        m_real += other.m_real;
        m_imag += other.m_imag;
        return *this;
    }

    ComplexNumber operator-(const ComplexNumber& other) const {
        return ComplexNumber(m_real - other.m_real, m_imag - other.m_imag);
    }

    ComplexNumber& operator-=(const ComplexNumber& other) {
        m_real -= other.m_real;
        m_imag -= other.m_imag;
        return *this;
    }

    ComplexNumber operator+(double value) const {
        return ComplexNumber(m_real + value, m_imag);
    }

    ComplexNumber& operator+=(double value) {
        m_real += value;
        return *this;
    }

    ComplexNumber operator-(double value) const {
        return ComplexNumber(m_real - value, m_imag);
    }

    ComplexNumber& operator-=(double value) {
        m_real -= value;
        return *this;
    }

    ComplexNumber operator*(const ComplexNumber& other) const {
        double tempReal = m_real * other.m_real - m_imag * other.m_imag;
        double tempImag = m_real * other.m_imag + other.m_real * m_imag;
        return ComplexNumber(tempReal, tempImag);
    }

    ComplexNumber& operator*=(const ComplexNumber& other) {
        double tempReal = m_real * other.m_real - m_imag * other.m_imag;
        double tempImag = m_real * other.m_imag + other.m_real * m_imag;
        m_real = tempReal;
        m_imag = tempImag;
        return *this;
}

    ComplexNumber operator/(const ComplexNumber& other) const {
        double denominator = (pow(other.m_real,2) + pow(other.m_imag,2));
        double tempReal = (m_real*other.m_real + m_imag*other.m_imag)/denominator;
        double tempImag = (-m_real*other.m_imag + other.m_real*m_imag)/denominator;
        return ComplexNumber(tempReal,tempImag);   
    }
    
    ComplexNumber& operator/=(const ComplexNumber& other) {
        double denominator = (pow(other.m_real,2) + pow(other.m_imag,2));
        double tempReal = (m_real*other.m_real + m_imag*other.m_imag)/denominator;
        double tempImag =  (-m_real*other.m_imag + other.m_real*m_imag)/denominator;
        m_real = tempReal;
        m_imag = tempImag;
        return *this;
    }

    bool operator==(const ComplexNumber& other) const {
        return (m_real == other.m_real) && (m_imag == other.m_imag);
    }

    bool operator!=(const ComplexNumber& other) const {
        return !(*this==other);
    }

    ComplexNumber& operator=(const ComplexNumber& other) {
    if (this != &other) { 
        m_real = other.m_real;
        m_imag = other.m_imag;
    }

    return *this;
}
    
    friend std::ostream& operator<<(std::ostream& os, const ComplexNumber& c) {
        os << c.m_real 
           << (c.m_imag >= 0 ? " + " : " - ") 
           << std::abs(c.m_imag) << "i";
        return os;
    }

};

