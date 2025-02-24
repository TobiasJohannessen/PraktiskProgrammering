#include "vec.h"

int main(){

    vec v1 = vec(1,2,3);
    vec v2 = vec(3,2,1);
    v1.print("v1 = ");
    v2.print("v2 = ");
    vec v3;
    v3 = v1 + v2;
    v3.print("v3 = v1 + v2 = ");
	v1 += v2;
    v1.print("v1 += v2, v1 = ");
    v1.set(1,2,3);
    return 0;
}