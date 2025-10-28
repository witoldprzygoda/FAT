#include "PPip_ID.h"
#include "src/manager.h"
#include "src/progressbar.h"
#include "data.h"
#include <TChain.h>
#include <TClass.h>
#include <TVector3.h>
#include <TDataMember.h> // Include this header
#include <iostream>
#include <iomanip>

ClassImp(PPip_ID);

PPip_ID::PPip_ID(const char *treeName) : Base_ID(nullptr, "PPip_ID")
{
    TChain *chain = new TChain(treeName, "");

    // Add here more files to the chain if needed
    //chain->Add("/lustre/hades/user/przygoda/pat3/out/exp/060/01/lep060_01_new.root");
    #include "h68_10.list"

    fTree = chain;
    PopulateLeafMap();
    SetBranchAddresses();
    Init();
}

// This method is a bit redundant, but it guarantees that ntuple structure is defined and stable already here.
// The ntuple is booked after the first fill() methond call.
// If you think you don't need it, you can remove it and fill ntuple straight in ProcessEntries() method.
// Init() is called in constructor
void PPip_ID::Init()
{
    HNtuple &r = *nt; // this is only to hande below by name, not pointer

    r["p_p"] = -1.0;
    r["pip_p"] = -1.0;
    r["ppip_m"] = -1.0;
    r["ppip_oa"] = -1.0;
    r["isBest"] = -1.0;

    r.fill(); // now the structure is frozen
}

// This is the main battle field where you should implement your analysis
void PPip_ID::ProcessEntries()
{
    HNtuple &r = *nt; // this is only to hande below by name, not pointer

    Long64_t nentries = fTree->GetEntries();
    // Limited to 10 for testing, uncomment GetEntries() above to process all entries

    for (Long64_t i = 0; i < nentries; ++i)
    {
	progressbar(i, nentries);
        fTree->GetEntry(i);
        //std::cout << "Entry " << i << ": p_p = " << p_p << std::endl;
        //std::cout << "Entry " << i << ": pip_p = " << pip_p << std::endl;

        // create four-vectors, remember to use pip_p_corr_pip, and so on

        TVector3 v0, v1, v2, v3;
        //v0.SetXYZ(p_p * sin(D2R * p_theta) * cos(D2R * p_phi), p_p * sin(D2R * p_theta) * sin(D2R * p_phi), p_p * cos(D2R * p_theta));
        //v1.SetXYZ(pip_p * sin(D2R * pip_theta) * cos(D2R * pip_phi), pip_p * sin(D2R * pip_theta) * sin(D2R * pip_phi), pip_p * cos(D2R * pip_theta));
        v0.SetXYZ(p_p_corr_p * sin(D2R * p_theta) * cos(D2R * p_phi), p_p_corr_p * sin(D2R * p_theta) * sin(D2R * p_phi), p_p_corr_p * cos(D2R * p_theta));
        v1.SetXYZ(pip_p_corr_pip * sin(D2R * pip_theta) * cos(D2R * pip_phi), pip_p_corr_pip * sin(D2R * pip_theta) * sin(D2R * pip_phi), pip_p_corr_pip * cos(D2R * pip_theta));

        p->SetVectM(v0, 938.27231);
        pip->SetVectM(v1, 139.56995);


	if (isBest==1) {
          eVertX->Fill(eVertReco_x);
          eVertY->Fill(eVertReco_y);
          eVertZ->Fill(eVertReco_z);
	}


        // Exmaple histograms filled
	//if (isBest==1 && eVertReco_z > -200. && eVertReco_z < 0.) {
	//if (isBest==1 && eVertReco_z > -140. && eVertReco_z < -100. && (((int)trigbit)&4096)) {
	if (isBest==1 && eVertReco_z > -130. && eVertReco_z < -110. && (((int)trigbit)&4096)) {
           //hp_p->Fill(p_p);
           //hp_theta->Fill(p_theta);
           //hp_phi->Fill(p_phi);
           //hp_theta_phi->Fill(p_theta, p_phi);
	   //banana_p->Fill(p_beta_new,p_p);
	   //banana_pip->Fill(pip_beta_new,pip_p);

      *p_LAB = *p;
      *pip_LAB = *pip;
      *n = *beam - *p - *pip;
      *n_LAB = *n;
      *deltaPP = *p + *pip;
      *deltaPP_LAB = *p + *pip;
      *deltaP = *beam - *p;
      *deltaP_LAB = *beam - *p;

      // pwa extension
      *p_pip = *p + *pip;
      *n_pip = *n + *pip;
      *pn =  *p + *n;

      // helicity frames
      *n_PPIP = *n;
      *p_NPIP = *p;
      *pip_PN = *pip;

      *pip_PPIP = *pip;
      *pip_NPIP = *pip;
      *n_PN = *n;

      *proj_PPIP = *proj;
      *proj_NPIP = *proj;
      *proj_PN = *proj;


      p->Boost(0.0, 0.0, -(*beam).Beta() );
      n->Boost(0.0, 0.0, -(*beam).Beta() );
      pip->Boost(0.0, 0.0, -(*beam).Beta() );
      deltaP->Boost(0.0, 0.0, -(*beam).Beta() );
      deltaPP->Boost(0.0, 0.0, -(*beam).Beta() );

       // boosts
       pip_PPIP->Boost( -(*p_pip).BoostVector() );
       n_PPIP->Boost( -(*p_pip).BoostVector() );
       pip_NPIP->Boost( -(*n_pip).BoostVector() );
       p_NPIP->Boost( -(*n_pip).BoostVector() );
       n_PN->Boost( -(*pn).BoostVector() );
       pip_PN->Boost( -(*pn).BoostVector() );
       // projectile in PPIP and NPIP and PN frames
       proj_PPIP->Boost( -(*p_pip).BoostVector() );
       proj_NPIP->Boost( -(*n_pip).BoostVector() );
       proj_PN->Boost( -(*pn).BoostVector() );

       double m2_inv_deltaPP = deltaPP->M2();
       double m_inv_deltaPP = deltaPP->M() / 1000.;
       double deltaPP_CM_cosTheta = deltaPP->CosTheta();

      double ang_cut = 9.;

      double close_cut = 5.;
      double nonfit_close_cut = -5.;

      NoLeptonPip = !( pip_oa_lept < close_cut && pip_oa_lept > nonfit_close_cut );
      NoHadronPip = !( pip_oa_hadr < close_cut && pip_oa_hadr > nonfit_close_cut );
      NoLeptonP  = !( p_oa_lept  < close_cut && p_oa_lept  > nonfit_close_cut );
      NoHadronP  = !( p_oa_hadr  < close_cut && p_oa_hadr  > nonfit_close_cut );
      NoHadronPip = 1;
      NoHadronPip = 1;
      NoLeptonP  = 1;
      NoHadronP = 1;


            //PiPlusProton = NoLeptonPip && NoHadronPip  &&  NoLeptonP && NoHadronP && insideTarget;
            PiPlusProton = NoLeptonPip && NoHadronPip  &&  NoLeptonP && NoHadronP;


           r["p_p"] = p_p;
           r["pip_p"] = pip_p;
           r["ppip_m"] = p_pip->M();
           r["ppip_oa"] = R2D * Manager::openingangle(*p, *pip);
           r["isBest"] = isBest;

           r.fill();



	    double EFF = 1, WEIGHT = 1, Qfactor = 1, Afactor = 1;

       if ( m_inv_deltaPP >= 1.0 && m_inv_deltaPP <= 1.8)
       {
            //int masa_id = static_cast< int >( ( 2.*(m_inv_deltaPP  - 1.0) )*25. );
            int masa_id = static_cast< int >( (deltaPP->M()/1000. - 1.0) * 25 / 0.8 );
            int cos_id  = static_cast< int >( (deltaPP_CM_cosTheta + 1.)*10. );
            Qfactor = Qfactor_tab[ masa_id  ] [ cos_id  ];
       }
       Qfactor = 1;

	         //if (PiPlusProton && (((int)s.trigbit)&32) && s.trigdownscaleflag==1 &&  s.isBest>0
	         if (true /*PiPlusProton */
      //if (PiPlusProton && (((int)s.trigbit)&16) && s.trigdownscaleflag==1 &&  s.isBest>0
      //if (PiPlusProton && (((int)s.trigbit)&8) && s.trigdownscaleflag==1 &&  s.isBest>0
//          && s.pip_z < -10. && s.pip_z > -90. && s.pip_p > 50. && s.p_z < -10. && s.p_z > -90. && s.p_p > 50.
          /* EXPLICIT ASK for the trigger M2-SMART_OPPOSITE */
         // && TMath::Abs( s.p_sector - s.pip_sector ) == 3 && ( s.p_system==0 || s.pip_system==0 )
          //&& TMath::Abs( p_sector_calc - pip_sector_calc ) == 3 && ( p_system_calc==0 || pip_system_calc==0 )
          //&& WEIGHT > 0.0
         )
      {

         double m_n = n->M() / 1000.;
         double m_n2 = n->M2() / 1000000.;
         double m_p = p->M() / 1000.;
         double m_p2 = p->M2() / 1000000.;
         double m_pip = pip->M() / 1000.;
         double m_pip2 = pip->M2() / 1000000.;
	 //std::cout << " m_n " << m_n << std::endl;
         //if ( m_n > 0.878 && m_n < 0.998 ) // DIRECT WINDOW (0.938+-0.060)a- DEFAULT condition
         if ( m_n > 0.899 && m_n < 0.986 ) // fit to all data
         {

            // 4-momentum transfer
            // t1 p - D++ forward
            // t1 targ - D++ backward
            // t2 p - D+ forward
            // t2 targ - D+ backward
	    
	                 //cout << " WEIGHT " << WEIGHT << endl;
             exp_sig->Fill ( deltaPP->M()/1000., deltaPP->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );

             // A
             pwa_pip_costh->Fill( pip->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );
             pwa_p_costh->Fill( p->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );
             pwa_n_costh->Fill( n->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );

             // B
             pwa_pip_p->Fill( pip_LAB->P()/1000., EFF*WEIGHT*Qfactor*Afactor );
             pwa_p_p->Fill( p_LAB->P()/1000., EFF*WEIGHT*Qfactor*Afactor );
             pwa_n_p->Fill( n_LAB->P()/1000., EFF*WEIGHT*Qfactor*Afactor );

             // C
             pwa_ppip_m->Fill( p_pip->M()/1000., EFF*WEIGHT*Qfactor*Afactor );
             pwa_npip_m->Fill( n_pip->M()/1000., EFF*WEIGHT*Qfactor*Afactor );
             pwa_pn_m->Fill( pn->M()/1000., EFF*WEIGHT*Qfactor*Afactor );

             // E
             pwa_pip_helicity->Fill( TMath::Cos( Manager::openingangle( *pip_PPIP, *n_PPIP ) ), EFF*WEIGHT*Qfactor*Afactor );
             pwa_pipn_helicity->Fill( TMath::Cos( Manager::openingangle( *pip_NPIP, *p_NPIP ) ), EFF*WEIGHT*Qfactor*Afactor );
             pwa_n_helicity->Fill( TMath::Cos( Manager::openingangle( *n_PN, *pip_PN ) ), EFF*WEIGHT*Qfactor*Afactor );

             // F

             // bad GJ
             //pwa_pip_gj->Fill( pip_PPIP->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );
             //pwa_pipn_gj->Fill( pip_NPIP->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );
             //pwa_n_gj->Fill( n_PN->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );
             // with projectile
             pwa_pip_gj->Fill( TMath::Cos( Manager::openingangle( *pip_PPIP, *proj_PPIP ) ), EFF*WEIGHT*Qfactor*Afactor );
             pwa_pipn_gj->Fill( TMath::Cos( Manager::openingangle( *pip_NPIP, *proj_NPIP ) ), EFF*WEIGHT*Qfactor*Afactor );
             pwa_n_gj->Fill( TMath::Cos( Manager::openingangle( *n_PN, *proj_PN ) ), EFF*WEIGHT*Qfactor*Afactor );


            // D++
            mass_deltaPP->Fill( deltaPP->M() / 1000. , EFF*WEIGHT*Qfactor*Afactor );
            cos_theta_deltaPP->Fill( deltaPP->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );

            // D+
            mass_deltaP->Fill( deltaP->M() / 1000. , EFF*WEIGHT*Qfactor*Afactor );
            cos_theta_deltaP->Fill( deltaP->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );

            // p
            cos_theta_p->Fill( p->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );

            // n
            cos_theta_n->Fill( n->CosTheta(), EFF*WEIGHT*Qfactor*Afactor );


	 }

	 mass_n->Fill( m_n , EFF*WEIGHT*Qfactor*Afactor );
         mass_p->Fill( m_p , EFF*WEIGHT*Qfactor*Afactor );
         mass_pip->Fill( m_pip , EFF*WEIGHT*Qfactor*Afactor );
         int costh_id = 0, mass_id = 0;
         // D++ cos theta
         costh_id = static_cast< int >( (deltaPP->CosTheta() + 1.)*20. );
         //mass_n_tab[ costh_id ]->Fill( m_n , EFF*WEIGHT );
         // D++ inv mass
         if ( deltaPP->M() / 1000. > 0.8 && deltaPP->M() / 1000. < 1.8 )
         {
            mass_id = static_cast< int >( ( deltaPP->M() / 1000. - 0.8 )*100. );
            // cout << " --- DPP masa " << deltaPP->M() / 1000. << " indeks = " << mass_id << endl;
            mass_p_tab[ mass_id ]->Fill( m_n , EFF*WEIGHT );
         }
         // D+ cos theta
         costh_id = static_cast< int >( (deltaP->CosTheta() + 1.)*20. );
         //mass_nn_tab[ costh_id ]->Fill( m_n , EFF*WEIGHT );
         // D+ inv mass
         if ( deltaP->M() / 1000. > 0.8 && deltaP->M() / 1000. < 1.8 )
         {
            mass_id = static_cast< int >( ( deltaP->M() / 1000. - 0.8 )*100. );
            // cout << "DP masa " << deltaPP->M() / 1000. << " indeks = " << mass_id << endl;
            mass_pp_tab[ mass_id ]->Fill( m_n , EFF*WEIGHT );
         }

	 // 2-dim
         int masa_id, cos_id;
         if ( deltaPP->M()/1000. >= 1.0 && deltaPP->M()/1000. <= 1.8 ) {
	    masa_id = static_cast< int >((deltaPP->M()/1000. - 1.0) * 25 / 0.8);
            cos_id  = static_cast< int >( (deltaPP->CosTheta() + 1.)*10. );
            miss_mass_tab[masa_id][cos_id] -> Fill ( m_n , EFF*WEIGHT );
         }

      } // if true




	}
    }
}

// Do not modify this method
void PPip_ID::SetBranchAddresses()
{
    TClass *cls = TClass::GetClass(typeid(*this));
    TIter next(cls->GetListOfDataMembers());
    TDataMember *member;
    while ((member = (TDataMember *)next()))
    {
        const char *memberName = member->GetName();
        TLeaf *leaf = leafMap[memberName];
        if (leaf)
        {
            void *memberAddr = (char *)this + member->GetOffset();
            fTree->SetBranchAddress(memberName, memberAddr);
        }
    }
}
