#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <complex>
#include <cmath>

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

cpp_dec_float_50* Sk(float* Ejs , int q , std::complex<float> D, int kmax){
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
            p_coeffs[j] += -1.0*Ejs[i]*p_coeffs[j+1];
        }
    }
    // Adding the (-1)^q D to the original polynomial:
    p_coeffs[0] = p_coeffs[0] + pow(-1.0 , q)*real(D);

    //Computing the derivative of p(x) (i.e. p'(E)):
    cpp_dec_float_50 p_der_coeffs[q];
    for(int i = 0; i < q; i++){
        p_der_coeffs[i] = (i+1)*p_coeffs[i+1];
    }

    //Doing long division of p'(x)/p(x) to find the coefficients S_k
    //Recall that: p'(x)/p(x) = \sum_{k=0}^\infty \frac{S_k}{x^{k+1}}
    cpp_dec_float_50* S_k = new cpp_dec_float_50[kmax]; // S_k is a pointer to an array of multiprecision floats
    cpp_dec_float_50 divided[q] , der_array[q+1];
    cpp_dec_float_50 a;
    for(int i=0; i < q ; i++){
        divided[i] = p_der_coeffs[i];
    }
    for(int i=0; i < kmax; i++){
        a = divided[q-1]/p_coeffs[q];
        for(int j=0; j < q+1; j++){
            if(j == 0){
                der_array[j] = -1.0*a*p_coeffs[j];
            } else{
                der_array[j] = -1.0*a*p_coeffs[j] + divided[j-1];
            }
        }
        S_k[i] = a;
        for(int j = 0; j < q; j++){
            divided[j] = der_array[j];
        }
    }
    return S_k;
}

int main()
{
    // Testing the polynomial coefficient generator:
    float Energies[] = {-1.67, 2.5, 4.87, 7.8};
    int q = sizeof(Energies)/sizeof(float);
    int kmax = 6;
    // ---------------------------------------------
    std::complex<float> D = {1.5, 0};
    cpp_dec_float_50* S_k_array = Sk(Energies , q , D , kmax);
    for(int i=0; i<kmax; i++){
        std::cout<<"S_k["<<i<<"] = ";
        std::cout<< S_k_array[i] <<std::endl;
    }
    delete[] S_k_array;
    return 0;
}