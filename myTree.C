#define myTree_cxx

#include "myTree.h"

void myTree::makeFlatTree() {
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  cout <<"[INFO] making the flat tree for the unfolding"<<endl;

  string trUnfFileName = Form("tree_%s_PP_jetR%d_AccEff_JEC.root",isMC?Form("MC%dS",nS):"DATA",jetR04?4:3);

  TFile* trUnfFile = new TFile (trUnfFileName.c_str(),"RECREATE");
  TTree* trUnf = new TTree ("treeForUnfolding","tree used for the unfolding");

  Int_t evtNb; Int_t centr; Float_t jp_pt; Float_t jp_rap; Float_t jp_eta; Float_t jp_mass; Float_t jp_phi; Float_t jp_l; 
  Float_t jp_gen_pt; Float_t jp_gen_rap; Float_t jp_gen_eta; Float_t jp_gen_phi; 
  Float_t jt_pt; Float_t jt_pt_noZJEC; Float_t jt_pt_jettyCorr; Float_t jt_pt_JEU_Up; Float_t jt_pt_JEU_Down; Float_t jt_rap; Float_t jt_eta; Float_t jt_phi; 
  Float_t jt_ref_pt; Float_t jt_ref_rap; Float_t jt_ref_eta; Float_t jt_ref_phi; Float_t jt_pt_genZJEC;
  Float_t z; Float_t gen_z; Float_t z_jettyCorr;
  Float_t corr_AccEff; Float_t pt_hat; Float_t corr_ptw;

  cout<< "[INFO] Creating the tree to use in the unfolding" << endl;
  trUnf->Branch("evtNb", &evtNb, "evtNb/I");
  trUnf->Branch("centr", &centr, "centr/I");
  trUnf->Branch("jp_pt", &jp_pt, "jp_pt/F");
  trUnf->Branch("jp_rap", &jp_rap, "jp_rap/F");
  trUnf->Branch("jp_eta", &jp_eta, "jp_eta/F");
  trUnf->Branch("jp_phi", &jp_phi, "jp_phi/F");
  trUnf->Branch("jp_mass", &jp_mass, "jp_mass/F");
  trUnf->Branch("jp_l", &jp_l, "jp_l/F");
  trUnf->Branch("jt_pt", &jt_pt, "jt_pt/F");
  trUnf->Branch("jt_pt_noZJEC", &jt_pt_noZJEC,"jt_pt_noZJEC/F");
  trUnf->Branch("jt_pt_jettyCorr", &jt_pt_jettyCorr,"jt_pt_jettyCorr/F");
  trUnf->Branch("jt_pt_JEU_Up", &jt_pt_JEU_Up,"jt_pt_JEU_Up/F");
  trUnf->Branch("jt_pt_JEU_Down", &jt_pt_JEU_Down,"jt_pt_JEU_Down/F");
  trUnf->Branch("jt_rap", &jt_rap, "jt_rap/F");
  trUnf->Branch("jt_eta", &jt_eta, "jt_eta/F");
  trUnf->Branch("jt_phi", &jt_phi, "jt_phi/F");
  trUnf->Branch("z", &z, "z/F");
  trUnf->Branch("z_jettyCorr", &z_jettyCorr, "z_jettyCorr/F");
  trUnf->Branch("corr_AccEff", &corr_AccEff, "corr_AccEff/F");
  if (isMC) {
    trUnf->Branch("corr_ptw", &corr_ptw, "corr_ptw/F");
    trUnf->Branch("pt_hat", &pt_hat, "pt_hat/F");
    trUnf->Branch("jp_gen_pt", &jp_gen_pt, "jp_gen_pt/F");
    trUnf->Branch("jp_gen_rap", &jp_gen_rap, "jp_gen_rap/F");
    trUnf->Branch("jp_gen_eta", &jp_gen_eta, "jp_gen_eta/F");
    trUnf->Branch("jp_gen_phi", &jp_gen_phi, "jp_gen_phi/F");
    trUnf->Branch("jt_ref_pt", &jt_ref_pt, "jt_ref_pt/F");
    trUnf->Branch("jt_ref_rap", &jt_ref_rap, "jt_ref_rap/F");
    trUnf->Branch("jt_ref_eta", &jt_ref_eta, "jt_ref_eta/F");
    trUnf->Branch("jt_ref_phi", &jt_ref_phi, "jt_ref_phi/F");
    trUnf->Branch("jt_pt_genZJEC", &jt_pt_genZJEC, "jt_pt_genZJEC/F");
    trUnf->Branch("gen_z", &gen_z, "gen_z/F");
  }

  //get the JECs
  vector<string> jecFileName;
  string jeuFileName;  
  jecFileName.push_back(Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Input/JECDatabase/textFiles/Spring18_ppRef5TeV_V4_%s/Spring18_ppRef5TeV_V4_%s_L2Relative_AK%dPF.txt",isMC?"MC":"DATA",isMC?"MC":"DATA",jetR04?4:3));
  jeuFileName = Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Input/JECDatabase/textFiles/Spring18_ppRef5TeV_V4_%s/Spring18_ppRef5TeV_V4_%s_Uncertainty_AK%dPF.txt","DATA","DATA",jetR04?4:3);
    if (!isMC)
      jecFileName.push_back(Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Input/JECDatabase/textFiles/Spring18_ppRef5TeV_V4_%s/Spring18_ppRef5TeV_V4_%s_L2L3Residual_AK%dPF.txt",isMC?"MC":"DATA",isMC?"MC":"DATA",jetR04?4:3));
JetCorrector JEC(jecFileName);
JetUncertainty JEU(jeuFileName);

  // starting to process the trees
  cout << "[INFO] Starting to process " << nentries << " nentries" << endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (jentry%1000000==0) cout <<"[INFO] entry = "<<jentry<<"/"<<nentries<<endl;

      //global event quantities
      evtNb = evt;
      centr = 0;
      pt_hat = pthat;
      corr_ptw = Gen_weight;

      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) {
	//cout <<"in QQ loop:"<<iQQ<<"/"<<Reco_QQ_size<<endl;
	//reset all quantities
	z=-99;
	gen_z=-99;
	TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	//cout <<"RecoQQ4mom made"<<endl;

	//cout <<"Reco_QQ_mupl_idx[iRecoQQ] = "<<Reco_QQ_mupl_idx[iQQ]<<"/"<<Reco_mu_size<<endl;
	//cout <<"Reco_QQ_mumi_idx[iRecoQQ] = "<<Reco_QQ_mumi_idx[iQQ]<<"/"<<Reco_mu_size<<endl;

	//if (Reco_QQ_mupl_idx[iQQ] < 0 || Reco_QQ_mumi_idx[iQQ] < 0) continue; //temporary
	//apply the QQ filters
	if (!(pPAprimaryVertexFilter && pBeamScrapingFilter)) continue;
	//cout <<"passed global event selection"<<endl; 
	if (!areMuonsInAcceptance2018(iQQ)) continue;
	//cout <<"passed acceptance"<<endl;
	if (!passQualityCuts2018(iQQ)) continue;
	//cout <<"passed ID for triggerIndex_PP = "<<triggerIndex_PP<<endl;
	if (!isTriggerMatch(iQQ, triggerIndex_PP)) continue;
	//cout <<"passed trigger"<<endl;
	if (Reco_QQ_sign[iQQ]!=0) continue;
	if (isMC && !isMatchedDiMuon(iQQ)) continue;
	//cout <<"passed all the selection"<<endl;

	//fill the QQ quantities
	jp_pt = RecoQQ4mom->Pt();
	jp_rap = RecoQQ4mom->Rapidity();
	jp_eta = RecoQQ4mom->Eta();
	jp_phi = RecoQQ4mom->Phi();
	jp_mass = RecoQQ4mom->M();
	jp_l = Reco_QQ_ctau[iQQ];
	corr_AccEff = 1.0/getAccEff(RecoQQ4mom->Rapidity(),RecoQQ4mom->Pt());

	if (isMC) {
	  //cout <<"Reco_QQ_whichGen[iQQ] = "<<Reco_QQ_whichGen[iQQ]<<"/"<<Gen_QQ_size<<endl;
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(Reco_QQ_whichGen[iQQ]);
	  jp_gen_pt = GenQQ4mom->Pt();
	  jp_gen_rap = GenQQ4mom->Rapidity();
	  jp_gen_eta = GenQQ4mom->Eta();
	  jp_gen_phi = GenQQ4mom->Phi();
	}

	//look for a jet
	bool jetFound = false;
	double drmin= 0.5;

	for (int iJet=0; iJet<nref; iJet++) {
	  //cout <<"in jet loop"<<endl;
	  if (jetFound) continue;
	  JEC.SetJetPT(rawpt[iJet]);
	  JEC.SetJetEta(jteta[iJet]);
	  JEC.SetJetPhi(jtphi[iJet]);
	  JEC.SetRho(1);
	  JEC.SetJetArea(jtarea[iJet]);
	  jt_pt_noZJEC = JEC.GetCorrectedPT();
	  
	  JEU.SetJetPT(jt_pt_noZJEC);
	  JEU.SetJetEta(jteta[iJet]);
	  JEU.SetJetPhi(jtphi[iJet]);
	  
	  TLorentzVector v_jet;
	  v_jet.SetPtEtaPhiM(jt_pt_noZJEC, jteta[iJet], jtphi[iJet], jtm[iJet]);
	  if (RecoQQ4mom->DeltaR(v_jet)<=drmin) {
	    drmin = RecoQQ4mom->DeltaR (v_jet);
	    jt_pt = jecCorr(jt_pt_noZJEC, rawpt[iJet], RecoQQ4mom->Pt());
	    z = jp_pt/jt_pt;
	    if (z >= 1 && z <= 1.0001) z = 0.9999;
	    jt_pt_JEU_Down = jeuCorr(jt_pt, z, -1*(JEU.GetUncertainty().first));//jt_pt * (1 - JEU.GetUncertainty().first);
	    jt_pt_JEU_Up = jeuCorr(jt_pt, z, JEU.GetUncertainty().second);//jt_pt * (1 + JEU.GetUncertainty().second);
	    jetFound = true;
	  } // end of dR condition
	  jt_rap = jty[iJet];
	  jt_eta = jteta[iJet];
	  jt_phi = jtphi[iJet];

	  if (rawpt[iJet]>jp_pt) {
	    JEC.SetJetPT(rawpt[iJet]-jp_pt); 
	    jt_pt_jettyCorr = jp_pt+JEC.GetCorrectedPT();
	  }
	  else jt_pt_jettyCorr = jp_pt;
	  z_jettyCorr = jp_pt/jt_pt_jettyCorr;

	  if (isMC) {
	    jt_ref_pt = refpt[iJet];
	    jt_ref_rap = refy[iJet];
	    jt_ref_eta = refeta[iJet];
	    jt_ref_phi = refphi[iJet];
	    gen_z = jp_gen_pt/jt_ref_pt;
	    if (gen_z >= 1 && gen_z <= 1.0001) gen_z = 0.9999;
	    jt_pt_genZJEC = zjecCorr(jt_pt_noZJEC, rawpt[iJet], gen_z);
	  }
	} //end of iJet loop
	trUnf->Fill();
      } //end of iQQ loop
   } //end of jentry loop
   cout<<"[INFO] Saving the tree for the unfolding" << endl;
   trUnfFile->Write();
   //trUnfFile->cd();
   //trUnf->Write();
   trUnfFile->Close(); delete trUnfFile;
}
