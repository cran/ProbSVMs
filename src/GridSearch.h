#ifndef GridSearch_H
#define GridSearch_H

#include <vector>
#include <cmath>

template <class Evaluate, class RetClass>
double CSearch(const int Cgridinlev, const int nloops, const double powerbase, Evaluate& f,
  RetClass* cstresp=NULL, RetClass* bstresp=NULL, double C0=1.)
{
   std::vector< std::vector <double> > validres;
   validres.resize(Cgridinlev);
   for (int i=0;i<Cgridinlev;++i) validres[i].resize(nloops);

   std::vector<double> Cgrid(Cgridinlev);
   int b, minb(Cgridinlev/2), meadind(minb+1), bestb(meadind);
   double C,bestC(1.),bestres(-std::numeric_limits<double>::max()),curres;

   for (int loop=0;loop<nloops;++loop) {
     for (int bind1=0,bind2=-minb; bind1<Cgridinlev; ++bind1,++bind2)
     {
       if (loop==0) {
         if (nloops==1) Cgrid[bind1] = C0*pow(powerbase,bind2);
         else Cgrid[bind1] = C0*pow(powerbase,bind2*meadind);
       } else {
         if (loop==1) Cgrid[bind1] = pow(powerbase,bind2)*bestC;
         else if (loop==2) Cgrid[bind1] = pow(powerbase,bind2/meadind)*bestC;
       }
     }
     for (int bind1=0;bind1<Cgridinlev;++bind1) {
       if (bind1==0) b=meadind;
       else if (bind1<meadind) b=bind1;
            else b=bind1+1;
       if (b==meadind && loop>0) {
         validres[meadind-1][loop] = validres[bestb-1][loop-1];
         bestb = meadind;        

       } else {
        C = Cgrid[b-1];
        curres = validres[b-1][loop] = f(C,cstresp);        
// Rprintf("C = %f curres = %f previous best result = %f\n",C,curres,bestres);        
        if (curres > bestres) {
           bestb = b;
           bestC = C;
           bestres = curres;
           if (cstresp!=NULL) *bstresp = *cstresp;
         }

       }
     }

   }

   return bestC;
}

#endif

