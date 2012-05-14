#include <math.h>

template<unsigned M, unsigned N, unsigned B, unsigned A>
struct SinCosSeries {
   static double value() {
      return 1-(A*M_PI/B)*(A*M_PI/B)/M/(M+1)
               *SinCosSeries<M+2,N,B,A>::value();
   }
};
 
template<unsigned N, unsigned B, unsigned A>
struct SinCosSeries<N,N,B,A> {
   static double value() { return 1.; }
};
 
template<unsigned B, unsigned A, typename T=double>
struct Sin;
 
template<unsigned B, unsigned A>
struct Sin<B,A,float> {
   static float value() {
      return (A*M_PI/B)*SinCosSeries<2,24,B,A>::value();
   }
};
template<unsigned B, unsigned A>
struct Sin<B,A,double> {
   static double value() {
      return (A*M_PI/B)*SinCosSeries<2,34,B,A>::value();
   }
};
 
template<unsigned B, unsigned A, typename T=double>
struct Cos;
 
template<unsigned B, unsigned A>
struct Cos<B,A,float> {
   static float value() {
      return SinCosSeries<1,23,B,A>::value();
   }
};
template<unsigned B, unsigned A>
struct Cos<B,A,double> {
   static double value() {
      return SinCosSeries<1,33,B,A>::value();
   }
};

int main() {
  
}