#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <time.h>
#include <math.h>
//-----------------------------------------------------------------------------

#define a7 78125
#define a11 48828125
#define a13 1220703125
#define a15 452807053
#define a17 582758085
#define mask 2147483647
#define then

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ( (x) < (y) ? (y) : (x) )

const int FALSE = 0;
const int TRUE = 1;

/** Class for generating random integer numbers etc. Modified code from
* previous class of Daniele Pretolani. Three routines can be used and a sign routine.
* \author Lars Relund Nielsen.
* \version 0.5
*/
class Random
{
public:

    /** Return a seed using the clock value.  */
   int Clock_seed();

    /** Initialization of random routine 1 (use X17). */
   void Init_len(int seed);

    /** Returns an integer, in the interval {mid,...,mad}. */
   int Int_length (int mid,int mad);

    /** Initialization of random routine 2 (use X11). */
   void Init_w(int seed);

    /** Returns an integer, in the interval {l,...,r}. */
   int Int_weight (int l,int r);

    /** Initialization of random sign routine (use X7). */
   void Init_sign(int seed);

    /** Returns an zero/one value. */
   int Sign ();

    /** Initialization of random routine 3 (use X15). */
   void Init_num(int seed);

    /** Returns an integer, in the interval {mid,...,mad}. */
   int Int_number (int mid, int mad);

    /** Calculate \f$Pr(X=x)\f$ when \f$X\sim bi(n,p)\f$.
    Modification of the code of Joe Nellis (mrknowitall@mtcrossroads.org) at http://www.codeproject.com.
    \pre Do no checking, i.e. assume that \f$p\in [0,1]\f$ and that \f$x\in [0,n]\f$.
    */
    double BinomPdf(double n, double p, double x);

private:

   int rnd7(int X);
   int rnd11(int X);
   int rnd13(int X);
   int rnd15(int X);
   int rnd17(int X);

   int IX7;
   int IX11;
   int IX13;
   int IX15;
   int IX17;
   int bitmask;
};

//-----------------------------------------------------------------------------

#endif
