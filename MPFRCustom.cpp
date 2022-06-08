#include "MPFRCustom.h"

MPFRCustom::MPFRCustom() {printf("blabla2\n");
    mpfr_init2(value,53);
}

MPFRCustom::MPFRCustom(const char * x) {printf("blabla1\n");
    MPFRCustom();printf("blabla3\n");
    *(this) = x;
}

MPFRCustom::MPFRCustom(const double x) {
    MPFRCustom();
    *(this) = x;
}

MPFRCustom::MPFRCustom(const int x) {
    MPFRCustom();
    *(this) = x;
}

MPFRCustom::~MPFRCustom() {
    mpfr_clear(value);
}

void MPFRCustom::operator=(const MPFRCustom x) {
    mpfr_set(value,x.value,MPFR_RNDN);
}

void MPFRCustom::operator=(const char * x) {printf("blabla4: %s\n",x);
    mpfr_set_str(value,"0.1",10,MPFR_RNDN);printf("blabla5\n");
}

void MPFRCustom::operator=(const double x) {
    mpfr_set_d(value,x,MPFR_RNDN);
}

void MPFRCustom::operator=(const int x) {
    mpfr_set_si(value,(long)x,MPFR_RNDN);
}

MPFRCustom operator+(const MPFRCustom x, const MPFRCustom y) {
    MPFRCustom result;
    mpfr_add(result.value,x.value,y.value,MPFR_RNDN);
    return result;
}

MPFRCustom & MPFRCustom::operator+=(const MPFRCustom x) {
    (*this) = (*this) + x;
    return *this;
}

MPFRCustom operator+(const MPFRCustom x, const double y) {
    MPFRCustom result;
    mpfr_add_d(result.value,x.value,y,MPFR_RNDN);
    return result;
}

MPFRCustom & MPFRCustom::operator+=(const double x) {
    (*this) = (*this) + x;
    return *this;
}

MPFRCustom operator+(const double x, const MPFRCustom y) {
    return y+x;
}

MPFRCustom operator-(const MPFRCustom x, const MPFRCustom y) {
    MPFRCustom result;
    mpfr_sub(result.value,x.value,y.value,MPFR_RNDN);
    return result;
}

MPFRCustom & MPFRCustom::operator-=(const MPFRCustom x) {
    (*this) = (*this) + x;
    return *this;
}

MPFRCustom operator-(const MPFRCustom x, const double y) {
    MPFRCustom result;
    mpfr_sub_d(result.value,x.value,y,MPFR_RNDN);
    return result;
}

MPFRCustom & MPFRCustom::operator-=(const double x) {
    (*this) = (*this) + x;
    return *this;
}

MPFRCustom operator-(const double x, const MPFRCustom y) {
    MPFRCustom result;
    mpfr_si_sub(result.value,x,y.value,MPFR_RNDN);
    return result;
}

MPFRCustom & MPFRCustom::operator-() {
    mpfr_neg(value,value,MPFR_RNDN);
    return *this;
}

MPFRCustom operator*(const MPFRCustom x, const MPFRCustom y) {
    MPFRCustom result;
    mpfr_mul(result.value,x.value,y.value,MPFR_RNDN);
    return result;
}

MPFRCustom & MPFRCustom::operator*=(const MPFRCustom x) {
    (*this) = (*this) * x;
    return *this;
}

MPFRCustom operator*(const MPFRCustom x, const double y) {
    MPFRCustom result;
    mpfr_mul_d(result.value,x.value,y,MPFR_RNDN);
    return result;
}

MPFRCustom & MPFRCustom::operator*=(const double x) {
    (*this) = (*this) * x;
    return *this;
}

MPFRCustom operator*(const double x, const MPFRCustom y) {
    return y*x;
}

MPFRCustom operator/(const MPFRCustom x, const MPFRCustom y) {
    MPFRCustom result;
    mpfr_div(result.value,x.value,y.value,MPFR_RNDN);
    return result;
}

MPFRCustom & MPFRCustom::operator/=(const MPFRCustom x) {
    (*this) = (*this) / x;
    return *this;
}

MPFRCustom operator/(const MPFRCustom x, const double y) {
    MPFRCustom result;
    mpfr_div_d(result.value,x.value,y,MPFR_RNDN);
    return result;
}

MPFRCustom & MPFRCustom::operator/=(const double x) {
    (*this) = (*this) / x;
    return *this;
}

MPFRCustom operator/(const double x, const MPFRCustom y) {
    MPFRCustom result;
    mpfr_d_div(result.value,x,y.value,MPFR_RNDN);
    return result;
}

MPFRCustom operator/(const int x, const MPFRCustom y) {
    MPFRCustom result;
    mpfr_si_div(result.value,x,y.value,MPFR_RNDN);
    return result;
}

MPFRCustom abs(const MPFRCustom x) {
    MPFRCustom result;
    mpfr_abs(result.value,x.value,MPFR_RNDN);
    return result;
}

MPFRCustom sqrt(const MPFRCustom x) {
    MPFRCustom result;
    mpfr_sqrt(result.value,x.value,MPFR_RNDN);
    return result;
}

MPFRCustom pow(const MPFRCustom x, const int y) {
    MPFRCustom result;
    mpfr_pow_si(result.value,x.value,y,MPFR_RNDN);
    return result;
}
