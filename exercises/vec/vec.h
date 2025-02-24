#ifndef HAVE_VEC_H
#define HAVE_VEC_H
#include<iostream>
#include<string>
#include<cmath>

struct vec{
    double x,y,z;
    vec(double x,double y,double z): x(x),y(y),z(z){} // parametrized constructor
    vec():vec(0,0,0){} // default constructor
    vec(const vec&)=default; // copy constructor
    vec(vec&&)=default; // move constructor
    ~vec()=default; // destructor
    vec& operator=(const vec&)=default; // copy assignment
    vec& operator=(vec&&)=default; // move assignment
    vec& operator+=(const vec& v){x+=v.x; y+=v.y; z+=v.z; return *this;};
    vec& operator-=(const vec& v){x-=v.x; y-=v.y; z-=v.z; return *this;};
    vec& operator*=(double c){x*=c; y*=c; z*=c; return *this;};;
    vec& operator/=(double c){x/=c; y/=c; z/=c; return *this;};;
    void set(double a,double b,double c){x=a;y=b;z=c;}
    void print(std::string s="") const{
        std::cout << s << x << " " << y << " " << z << std::endl; // for debugging
    }; 
    friend std::ostream& operator<<(std::ostream&, const vec&);
};
vec operator-(const vec& v){return vec(-v.x, -v.y, -v.z);};
vec operator-(const vec& v1, const vec& v2){return vec(v1.x -v2.x, v1.y - v2.y, v1.z - v2.z);};
vec operator+(const vec& v1, const vec& v2){return vec(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);};
vec operator*(const vec& v, double c){return vec(v.x*c, v.y*c, v.z*c);};
vec operator*(double c, const vec& v){return v*c;};
vec operator/(const vec& v, double c){return vec(v.x/c, v.y/c, v.z/c);};
bool approx(double a, double b, double acc=1e-9, double eps=1e-9){
    double diff = abs(a-b);
    double size = std::max(abs(a), abs(b));
    if (diff <= acc){
        return true;
    };
    if (diff/size <= eps){
        return true;
    };
    return false;

};
#endif