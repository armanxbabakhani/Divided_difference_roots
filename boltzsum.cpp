#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <complex>

using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_int;

cpp_int factorial(int n)
{
    cpp_int result = 1;
    for (int i = 2; i <= n; ++i)
    {
        result *= i;
    }
    return result;
}

// We need a function with the following inputs: 1 - array of floats E_js of length q, 2 - A complex D , 3 - An integer kmax

cpp_dec_float_50 Sk(float* Ejs , int q , std::complex<float> D, int kmax){
    // Computing the polynomial p(E):
    cpp_dec_float_50 p_coeffs[q+1];
    p_coeffs[0] = 1.0;
    for(int i = 0; i < q; i++){
        // Multiplication by E:
        for(int j = 0; j < i + 1 ; j++){
            p_coeffs[i+1 - j] = p_coeffs[i-j];
        }
        p_coeffs[0] = 0.0;
        // Multiplication by the root E_i:
        for(int j = 0; j < i + 1; j++){
            p_coeffs[j] += Ejs[i]*p_coeffs[j+1];
        }
    }
    //Computing the derivative of p(x) (i.e. p'(E)):
    cpp_dec_float_50 p_der_coeffs[q];
    for(int i = 0; i < q; i++){
        p_der_coeffs[i] = (i+1)*p_coeffs[i+1];
    }

    for(int i = 0; i < q+1; i++){
        std::cout << "p(x)[" << i << "] = " << p_coeffs[i] << std::endl;
    }
    for(int i = 0; i < q; i++){
        std::cout << "p'(x)[" << i << "] = " << p_der_coeffs[i] << std::endl;
    }
    return 0;
}

int main()
{
    // Testing the polynomial coefficient generator:
    float Energies[] = {-1.67, 2.5, 4.87, 7.8};
    int q = sizeof(Energies)/sizeof(float), kmax = 6;
    std::cout << "q is: " << q << std::endl;
    std::complex<float> D = {1.5, 0};
    cpp_dec_float_50 a = Sk(Energies , q , D , kmax);
    return 0;
}