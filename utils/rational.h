/**
 * @author Malkovsky Nikolay, april 2016
 * 
 * Class provides basic arithmetic operations with rational numbers.
 * Rational numbers represented as nominator / denominator with
 * gcd(nom, denom) = 1 and denom > 0.
 */

#ifndef RATIONAL_H
#define RATIONAL_H

#define __EUCLID_
//#define __TRIVIAL_GCD_


class Rational {
public:
	long long nom;
	long long denom;

	Rational();
	Rational(long long a, long long b);

   	Rational(long long a);

   	Rational operator + ( Rational other) ;
   	 Rational operator + ( int other) ;
   	 Rational operator + ( long long other) ; 
   	 friend Rational operator + ( int other,  Rational r) ;
   	 friend Rational operator + ( long long other,  Rational r) ; 
	
	 Rational &operator += ( Rational other) ;
   	 Rational &operator += ( int other) ;
   	 Rational &operator += ( long long other) ; 

   	 Rational operator - ( Rational other) ;
   	 Rational operator - ( int other) ;
   	 Rational operator - ( long long other) ;
   	 friend Rational operator - ( long long other,  Rational r) ;
   	 friend Rational operator - ( int other,  Rational r) ;
    
   	 Rational &operator -= ( Rational other) ;
   	 Rational &operator -= ( int other) ;
   	 Rational &operator -= ( long long other) ;
   	 
   	 Rational operator * ( Rational other) ;
   	 Rational operator * ( int other) ;
   	 Rational operator * ( long long other) ;
	 friend Rational operator * ( long long other,  Rational r) ;
	 friend Rational operator * ( int other,  Rational r) ;

	 Rational &operator *= ( Rational other) ;
   	 Rational &operator *= ( int other) ;
   	 Rational &operator *= ( long long other) ;
	 
	 Rational operator / ( Rational other) ;
   	 Rational operator / ( int other) ;
   	 Rational operator / ( long long other) ;
   	 friend Rational operator / ( long long other,  Rational r) ;
   	 friend Rational operator / ( int other,  Rational r) ;
    
     Rational &operator /= ( Rational other) ;
   	 Rational &operator /= ( int other) ;
   	 Rational &operator /= ( long long other) ;
   	 
   	bool operator < ( Rational other) ;
   	bool operator < ( int other) ;
   	bool operator < ( long long other) ;
   	friend bool operator < ( int other,  Rational r) ;
	friend bool operator < ( long long other,  Rational r) ;

    bool operator > ( const Rational &other) ;
   	bool operator > ( const int &other) ;
   	bool operator > ( const long long &other) ;
   	friend bool operator > ( const int &other,  const Rational &r) ;
	friend bool operator > ( const long long &other,  const Rational &r) ;

	bool operator <= ( const Rational &other) ;
   	bool operator <= ( const int &other) ;
   	bool operator <= ( const long long &other) ;
   	friend bool operator <= ( const int &other,  const Rational &r) ;
	friend bool operator <= ( const long long &other,  const Rational &r) ;
	
	bool operator >= ( const Rational &other) ;
   	bool operator >= ( const int &other) ;
   	bool operator >= ( const long long &other) ;
   	friend bool operator >= ( const int &other,  const Rational &r) ;
	friend bool operator >= ( const long long &other,  const Rational &r) ;

	bool operator == (const Rational &other) ;
	bool operator == ( const int &other) ;
   	bool operator == ( const long long &other) ;
   	friend bool operator == ( const int &other,  const Rational &r) ;
	friend bool operator == ( const long long &other,  const Rational &r) ;

	bool operator != (const Rational &other) ;
	bool operator != ( const int &other) ;
   	bool operator != ( const long long &other) ;
   	friend bool operator != ( const int &other,  const Rational &r) ;
	friend bool operator != ( const long long &other,  const Rational &r) ;

	double toDouble();   	
   	
   	Rational operator - () ;
};
#endif