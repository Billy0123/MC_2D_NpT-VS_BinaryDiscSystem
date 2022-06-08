#ifndef MPFRCUSTOM_H
#define MPFRCUSTOM_H
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <mpfr.h>


class MPFRCustom
{
public:
    mpfr_t value;
    mpfr_prec_t prec=256;

    MPFRCustom();
    MPFRCustom(const char *);
    MPFRCustom(const double);
    MPFRCustom(const int);
    ~MPFRCustom();
    void operator=(const MPFRCustom);
    void operator=(const char *);
    void operator=(const double);
    void operator=(const int);
    friend MPFRCustom operator+(const MPFRCustom,const MPFRCustom);
    MPFRCustom & operator+=(const MPFRCustom);
    friend MPFRCustom operator+(const MPFRCustom,const double);
    MPFRCustom & operator+=(const double);
    friend MPFRCustom operator+(const double,const MPFRCustom);
    friend MPFRCustom operator-(const MPFRCustom,const MPFRCustom);
    MPFRCustom & operator-=(const MPFRCustom);
    friend MPFRCustom operator-(const MPFRCustom,const double);
    MPFRCustom & operator-=(const double);
    friend MPFRCustom operator-(const double,const MPFRCustom);
    MPFRCustom & operator-();
    friend MPFRCustom operator*(const MPFRCustom,const MPFRCustom);
    MPFRCustom & operator*=(const MPFRCustom);
    friend MPFRCustom operator*(const MPFRCustom,const double);
    MPFRCustom & operator*=(const double);
    friend MPFRCustom operator*(const double,const MPFRCustom);
    friend MPFRCustom operator/(const MPFRCustom,const MPFRCustom);
    MPFRCustom & operator/=(const MPFRCustom);
    friend MPFRCustom operator/(const MPFRCustom,const double);
    MPFRCustom & operator/=(const double);
    friend MPFRCustom operator/(const double,const MPFRCustom);
    friend MPFRCustom operator/(const int,const MPFRCustom);

    friend MPFRCustom abs(const MPFRCustom);
    friend MPFRCustom sqrt(const MPFRCustom);
    friend MPFRCustom pow(const MPFRCustom,const int);
};

#endif // MPFRCUSTOM_H
