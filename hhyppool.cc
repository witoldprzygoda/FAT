#include "hades.h"
#include "hevent.h"
#include "heventheader.h"
#include "hhyppool.h"
#include <algorithm>


// ****************************************************************************
ClassImp(HHypPool)

//-----------------------------------------------------------------------------
HHypPool::HHypPool(HOutputFile *ptr) : HPool(ptr), HHypDataPool(), ptrFH(new HParticlePlayerHandle) 
{
   reset();
}

//-----------------------------------------------------------------------------
HHypPool::~HHypPool() 
{ 
   reset();  
}

   
//-----------------------------------------------------------------------------
void HHypPool::reset()
{
   hypNum.clear();

   std::for_each( hypCand.begin(), hypCand.end(), DelHypCand() ); 
   hypCand.clear();
}


//-----------------------------------------------------------------------------
void HHypPool::loop(HParticlePool& ref)
{
   objIt = objectList.begin();
   while ( objIt != objectList.end() )
   {
      count( *objIt, ref );
      combine( *objIt, ref );
      for (hypIt = hypSet.begin(); hypIt != hypSet.end(); ++hypIt)
      {
         addHypCand( (*objIt)->getName(), *hypIt );
      }
      ++objIt;
   }
}

//-----------------------------------------------------------------------------
void HHypPool::fill() 
{
   HParticle *ptrP = 0;
   HHypCandidate *ptrH = 0;
   std::string key;
   //EParticle eid;
   
   //eventData.set("bbb_test", 4.0 );
   eventData.update();
   //eventData.set("bbb2_test", 5.0 );
   
   if (ptrFile != 0)
   {
       // internal loop over all hyp patterns
       for (objIt = objectList.begin(); objIt != objectList.end(); ++objIt)
       {
       	  // loop over all hyps stored
          hypIter = hypCand.begin();
          while(hypIter != hypCand.end())
          {
       	     if ( (*objIt)->getName() == hypIter->first )
       	     {
       	        ptrH = hypIter->second;
				numId.clear();
       	        for (int i=0; i<ptrH->getSize(); ++i)
				{
			       ++numId[ ptrH->getPart(i)->getEId() ];	
				}

       	        for (int i=0; i<ptrH->getSize(); ++i)
                {
                   ptrP = ptrH->getPart(i);
                   key  = convertId( convertId( ptrP->getId() ) );
				   if ( ptrH->getNum( ptrP->getEId() ) > 1 )
				   {
				      int k = numId[ ptrH->getPart(i)->getEId() ]--;
				      (*objIt)->setPrefix(key+static_cast<char>(k+48)+std::string("_"));
				   }
				   else
				   {
				      (*objIt)->setPrefix(key+"_");
				   }
                   (*objIt)->setSuffix(); 

                   for (NtuplePairIter it = ptrP->begin(); it != ptrP->end(); ++it)
                   {
                      (*objIt)->set( it->first, it->second );
                   }
                   for (NtuplePairIter it = ptrP->getCandidate()->begin(); it != ptrP->getCandidate()->end(); ++it)
                   {
                      if (ptrP->isName(it->first) == false)
                      {
                         (*objIt)->set( it->first, it->second );
                      }
                   }
                } 
               
               (*objIt)->setPrefix();
               (*objIt)->setSuffix();
               for (NtuplePairIter it = eventData.begin(); it != eventData.end(); ++it)
               {
                  (*objIt)->set( it->first, it->second );
               }
               
               (*objIt)->fill(); 
       	     }
       	     ++hypIter;
          }
       } // end of patterns
   } // if ptrFile
   
   //eventData.set("bbb3_test", 6.0 );
}


//-----------------------------------------------------------------------------
int HHypPool::count(HPattern *ptrC, HParticlePool& refPPool) 
{
   for( ptrC->partNumIter = ptrC->partComb.begin(); ptrC->partNumIter != ptrC->partComb.end(); ++ptrC->partNumIter )
   {
      // reset number of combinations for a given particle species
      ptrC->partNumIter->second = 0;
   }

   ptrC->number = 1;
   for( ptrC->partNumIter = ptrC->partNum.begin(); ptrC->partNumIter != ptrC->partNum.end(); ++ptrC->partNumIter )
   {
      if (refPPool.isPart( ptrC->partNumIter->first ) == true)
      {
         // ROOT changed convention !!! of what Binomial returns
         ptrC->partComb[ ptrC->partNumIter->first ] = ( refPPool.getNum( ptrC->partNumIter->first ) < ptrC->partNumIter->second ) ? 0 :
#if COMB_REPETITION == 1
              static_cast<int>( TMath::Binomial( refPPool.getNum( ptrC->partNumIter->first ), ptrC->partNumIter->second )*TMath::Factorial(ptrC->partNumIter->second) );
#else
              static_cast<int>( TMath::Binomial( refPPool.getNum( ptrC->partNumIter->first ), ptrC->partNumIter->second ) );
#endif
         ptrC->number *= ptrC->partComb[ ptrC->partNumIter->first ];
      }
      else
      {
         ptrC->number = 0;
      }
   }
return ptrC->number;
}


//-----------------------------------------------------------------------------
void HHypPool::combine(HPattern *ptrC, HParticlePool& refPPool)
{
   int i1, i2, i3, i4, repetitionNr, combinitaionNr;
   EParticle pid;

   hypSet.clear();
   count(ptrC, refPPool);
   hypSet.resize( ptrC->number );
   hypIt = hypSet.begin();
   while (hypIt != hypSet.end())
   {
      *hypIt = new HHypCandidate( ptrC );
      ++hypIt;
   }
   // here all HParticleCandidate objects from HParticlePool should be combined according to pattern from mapN

   i1 = i2 = i3 = i4 = 0; // this is to avoid warning about unused variables
   MultiParticleIter mpIter;
   MultiParticleIter mpBeginIter;
   for (ptrC->partNumIter = ptrC->partComb.begin(); ptrC->partNumIter != ptrC->partComb.end(); ++ptrC->partNumIter)
   // loop over all particle species
   {
      pid = ptrC->partNumIter->first;
      combinitaionNr = ptrC->partNumIter->second;
      if (combinitaionNr > 0) // for given particle: if number of combinations > 0
      {
            repetitionNr = ptrC->number / combinitaionNr;
            for (int k=0; k < combinitaionNr; ++k) // all combinations
            {
               switch (ptrC->partNum[ pid ])
               {
                   case 1: i1 = ptrFH->map1[ 10 + refPPool.getNum( pid ) ][ k ][0];
                           mpBeginIter = refPPool.lower_bound( pid );
                           mpIter = mpBeginIter;
                           while (i1--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           break;

#if MAX_PARTICLES_IN_COMB > 1
                   case 2: i1 = ptrFH->map2[ 20 + refPPool.getNum( pid ) ][ k ][0];
                           i2 = ptrFH->map2[ 20 + refPPool.getNum( pid ) ][ k ][1];
                           mpBeginIter = refPPool.lower_bound( pid );
                           mpIter = mpBeginIter;
                           while (i1--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                           {
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           }
                           mpIter = mpBeginIter;
                           while (i2--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           break;
#endif

#if MAX_PARTICLES_IN_COMB > 2
                   case 3: i1 = ptrFH->map3[ 30 + refPPool.getNum( pid ) ][ k ][0];
                           i2 = ptrFH->map3[ 30 + refPPool.getNum( pid ) ][ k ][1];
                           i3 = ptrFH->map3[ 30 + refPPool.getNum( pid ) ][ k ][2];
                           mpBeginIter = refPPool.lower_bound( pid );
                           mpIter = mpBeginIter;
                           while (i1--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           mpIter = mpBeginIter;
                           while (i2--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           mpIter = mpBeginIter;
                           while (i3--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           break;
#endif

#if MAX_PARTICLES_IN_COMB > 3
                   case 4: i1 = ptrFH->map4[ 40 + refPPool.getNum( pid ) ][ k ][0];
                           i2 = ptrFH->map4[ 40 + refPPool.getNum( pid ) ][ k ][1];
                           i3 = ptrFH->map4[ 40 + refPPool.getNum( pid ) ][ k ][2];
                           i4 = ptrFH->map4[ 40 + refPPool.getNum( pid ) ][ k ][3];
                           mpBeginIter = refPPool.lower_bound( pid );
                           mpIter = mpBeginIter;
                           while (i1--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           mpIter = mpBeginIter;
                           while (i2--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           mpIter = mpBeginIter;
                           while (i3--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           mpIter = mpBeginIter;
                           while (i4--) ++mpIter;
                           for (int j=0; j < repetitionNr; ++j) // number of repetitions of combinations
                              *hypSet [ k + j*combinitaionNr ] += mpIter->second;
                           break;
#endif

               }
            }
      }
   }

}


