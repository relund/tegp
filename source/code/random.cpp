#include "random.hpp"

int Random::rnd7 (int IX )
{
  IX = IX * a7;
  if (IX < 0)
  then IX = IX & mask;
  return (IX);
}; 
 
int Random::rnd11 (int IX )
{
  IX = IX * a11;
  if (IX < 0)
  then IX = IX & mask;
  return (IX);
}; 
 
int Random::rnd13 (int IX )
{
  IX = IX * a13;
  if (IX < 0)
  then IX = IX & mask;
  return (IX);
}; 
 
int Random::rnd15 (int IX )
{
  IX = IX * a15;
  if (IX < 0)
  then IX = IX & mask;
  return (IX);
}; 
 
int Random::rnd17 (int IX )
{
  IX = IX * a17;
  if (IX < 0)
  then IX = IX & mask;
  return (IX);
}; 

int Random::Clock_seed()
{
  time_t times;
  int seed;
  unsigned int st;
  times = time(0);
  st = (unsigned int) times / 2;
  seed = (int) st;
  if ( (seed%2) == 0 ) seed = seed + 1;
  return(seed);
}


void Random::Init_len(int seed)
{
  IX17 = rnd17(seed);
}

int Random::Int_length (int mid,int mad)
{ /* returns an integer, in the interval [mid..mad] **/
   int j, range;
   double f; 
     range = mad - mid + 1;
     IX17 = rnd17(IX17);
     f = (double) (IX17 * 0.4656613E-9);
     j = mid + (int) ( f * range );
     return(j);
 }

void  Random::Init_w(int seed)
{
  IX11 = rnd11(seed);
}

int Random::Int_weight (int l,int r )
{  /* generates an integer in [l..r]  **/
   double f;
   int j, range;
     if (r < l) 
     then return(0);
     else if (r==l) return(l);
     range = r-l+1;
     IX11 = rnd11(IX11);
     f = (double) (IX11 * 0.4656613E-9);
     j = l + (int) ( f * range );
     return(j);
 }


void Random::Init_sign(int seed)
{
  IX7 = rnd7(seed);
  bitmask = 1<<15;
};

int Random::Sign ()
{ /* returns a zero/one value **/
     IX7 = rnd7(IX7);
     if (IX7 & bitmask)
     then return(1);
     else return(0);
 }


void Random::Init_num(int seed)
{
  IX15 = rnd15(seed);
}

int Random::Int_number(int mid,int mad)
{ /* returns an integer, in the interval [mid..mad] **/
   int j, range;
   double f; 
     range = mad - mid + 1;
     IX15 = rnd17(IX15);
     f = (double) (IX15 * 0.4656613E-9);
     j = mid + (int) ( f * range );
     return(j);
 }

double Random::BinomPdf(double n, double p, double x)
{
    if(n==0) return 0.0;
    // initialize some variables
    double result = 1.0;
    double q = 1.0 - p;     // komplement to p
    int range = 0, np =0, nq = 0, nnumer = 0, ndenom = 0;
    
    // simple cases
    if(x == 0) return pow(q,n);
    if(x == n) return pow(p,x);

    // reorder the factorials to account for cancellations
    // in numerator and denominator.
    if(x < n-x) range = x;		// n-x cancels out
    else range = n-x;	// x cancels out
    np = x;
    nq = n-x;
    ndenom = range;
    nnumer = n;
    
    while(np > 0 || nq > 0 || ndenom > 0 || nnumer >(n-range)){
        // If the result is greater than one we want to divide by 
        // a denominator digit or multiply by percentage p or q.
        // If we are out of numerator digits then finish multiplying
        // with our powers of p or q or dividing by a denom digit.
        if(result >= 1.0 || nnumer ==(n-range)){
            if(ndenom > 0){
                //m_resut *= (1.0/ndenom);
                result /= ndenom;
                --ndenom;
            }
            else if(nq > 0){
                result *= q;
                --nq;
            }
            else if(np > 0){
                result *= p;
                --np;
            }
            else {
                throw("Binomial Probability computation error- check success percentage between 0 and 1");
            } 
        }
        // If the result is less than one then we want to multiply
        // by a numerator digit. If we are out of denominator digits,
        // powers of p or powers of q then multiply rest of result 
        // by numerator digits.
        else if(result < 1.0 || np==0  ){
            if(nnumer >(n-range)){
                result *= nnumer;
                --nnumer;
            }
            else{
                throw("Binomial Probability computation error- unknown error");
            }
        }
        else{
            throw("Binomial Probability computation error- possible value infinity or NaN");
        }
    }
    return result; 
}