//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 19 13:48:11 2021 by ROOT version 6.14/04
// from TTree myTree/My TTree of dimuons
// found on file: /data_CMS/cms/mnguyen/upsJet/mc/merged_HiForestAOD.root
//////////////////////////////////////////////////////////

#ifndef myTree_h
#define myTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "../JpsiInJetsPbPb/Fitter/Macros/Utilities/JetCorrector.h"
#include "../JpsiInJetsPbPb/Fitter/Macros/Utilities/JetUncertainty.h"

using namespace std;

// Header file for the classes stored in the TTree if any.

class myTree {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.
  bool            isMC;
  bool            jetR04;
  int             nS;
  int             triggerIndex_PP = 3;

  // Declaration of leaf types: onia
  Int_t           Reco_QQ_size;
  Int_t           Reco_QQ_type[99];   //[Reco_QQ_size]
  Int_t           Reco_QQ_sign[99];   //[Reco_QQ_size]
  TClonesArray    *Reco_QQ_4mom;
  Int_t           Reco_QQ_mupl_idx[99];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_idx[99];   //[Reco_QQ_size]
  ULong64_t       Reco_QQ_trig[99];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_isCowboy[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctau[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctauErr[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_cosAlpha[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctau3D[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctauErr3D[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_cosAlpha3D[99];   //[Reco_QQ_size]
  Int_t           Reco_QQ_whichGen[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_dca[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_MassErr[99];   //[Reco_QQ_size]
  TClonesArray    *Reco_QQ_vtx;
  Int_t           Reco_QQ_Ntrk[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dxy_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxy_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dxyErr_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxyErr_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dz_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dz_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dzErr_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dzErr_muonlessVtx[99];   //[Reco_QQ_size]
  Int_t           Reco_mu_size;
  Int_t           Reco_mu_type[99];   //[Reco_mu_size]
  Int_t           Reco_mu_whichGen[99];   //[Reco_mu_size]
  Int_t           Reco_mu_SelectionType[99];   //[Reco_mu_size]
  Int_t           Reco_mu_charge[99];   //[Reco_mu_size]
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_mu_trig[99];   //[Reco_mu_size]
  Bool_t          Reco_mu_highPurity[99];   //[Reco_mu_size]
  Bool_t          Reco_mu_TrkMuArb[99];   //[Reco_mu_size]
  Bool_t          Reco_mu_TMOneStaTight[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nPixValHits[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nMuValHits[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nTrkHits[99];   //[Reco_mu_size]
  Float_t         Reco_mu_normChi2_inner[99];   //[Reco_mu_size]
  Float_t         Reco_mu_normChi2_global[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nPixWMea[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nTrkWMea[99];   //[Reco_mu_size]
  Int_t           Reco_mu_StationsMatched[99];   //[Reco_mu_size]
  Float_t         Reco_mu_dxy[99];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[99];   //[Reco_mu_size]
  Float_t         Reco_mu_dz[99];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[99];   //[Reco_mu_size]
  Float_t         Reco_mu_pt_inner[99];   //[Reco_mu_size]
  Float_t         Reco_mu_pt_global[99];   //[Reco_mu_size]
  Float_t         Reco_mu_ptErr_inner[99];   //[Reco_mu_size]
  Float_t         Reco_mu_ptErr_global[99];   //[Reco_mu_size]
  Float_t         Gen_weight;
  Float_t         Gen_pthat;
  Int_t           Gen_QQ_size;
  Int_t           Gen_QQ_type[99];   //[Gen_QQ_size]
  TClonesArray    *Gen_QQ_4mom;
  Int_t           Gen_QQ_momId[99];   //[Gen_QQ_size]
  Float_t         Gen_QQ_ctau[99];   //[Gen_QQ_size]
  Float_t         Gen_QQ_ctau3D[99];   //[Gen_QQ_size]
  Int_t           Gen_QQ_mupl_idx[99];   //[Gen_QQ_size]
  Int_t           Gen_QQ_mumi_idx[99];   //[Gen_QQ_size]
  Int_t           Gen_QQ_whichRec[99];   //[Gen_QQ_size]
  Int_t           Gen_mu_size;
  Int_t           Gen_mu_type[99];   //[Gen_mu_size]
  Int_t           Gen_mu_charge[99];   //[Gen_mu_size]
  TClonesArray    *Gen_mu_4mom;
  Int_t           Gen_mu_whichRec[99];   //[Gen_mu_size]

  // Declaration of leaf types: jet
  Int_t           evt;
  Float_t         b;
  Float_t         vx;
  Float_t         vy;
  Float_t         vz;
  Int_t           nref;
  Float_t         rawpt[99];   //[nref]
  Float_t         jtpt[99];   //[nref]
  Float_t         jteta[99];   //[nref]
  Float_t         jty[99];   //[nref]
  Float_t         jtphi[99];   //[nref]
  Float_t         jtpu[99];   //[nref]
  Float_t         jtm[99];   //[nref]
  Float_t         jtarea[99];   //[nref]
  Float_t         jtPfCHF[99];   //[nref]
  Float_t         jtPfNHF[99];   //[nref]
  Float_t         jtPfCEF[99];   //[nref]
  Float_t         jtPfNEF[99];   //[nref]
  Float_t         jtPfMUF[99];   //[nref]
  Int_t           jtPfCHM[99];   //[nref]
  Int_t           jtPfNHM[99];   //[nref]
  Int_t           jtPfCEM[99];   //[nref]
  Int_t           jtPfNEM[99];   //[nref]
  Int_t           jtPfMUM[99];   //[nref]
  Float_t         jttau1[99];   //[nref]
  Float_t         jttau2[99];   //[nref]
  Float_t         jttau3[99];   //[nref]
  Float_t         discr_jetID_cuts[99];   //[nref]
  Float_t         discr_jetID_bdt[99];   //[nref]
  Float_t         discr_fr01[99];   //[nref]
  Float_t         trackMax[99];   //[nref]
  Float_t         trackSum[99];   //[nref]
  Int_t           trackN[99];   //[nref]
  Float_t         trackHardSum[99];   //[nref]
  Int_t           trackHardN[99];   //[nref]
  Float_t         chargedMax[99];   //[nref]
  Float_t         chargedSum[99];   //[nref]
  Int_t           chargedN[99];   //[nref]
  Float_t         chargedHardSum[99];   //[nref]
  Int_t           chargedHardN[99];   //[nref]
  Float_t         photonMax[99];   //[nref]
  Float_t         photonSum[99];   //[nref]
  Int_t           photonN[99];   //[nref]
  Float_t         photonHardSum[99];   //[nref]
  Int_t           photonHardN[99];   //[nref]
  Float_t         neutralMax[99];   //[nref]
  Float_t         neutralSum[99];   //[nref]
  Int_t           neutralN[99];   //[nref]
  Float_t         hcalSum[99];   //[nref]
  Float_t         ecalSum[99];   //[nref]
  Float_t         eMax[99];   //[nref]
  Float_t         eSum[99];   //[nref]
  Int_t           eN[99];   //[nref]
  Float_t         muMax[99];   //[nref]
  Float_t         muSum[99];   //[nref]
  Int_t           muN[99];   //[nref]
  Int_t           beamId1;
  Int_t           beamId2;
  Float_t         pthat;
  Float_t         refpt[99];   //[nref]
  Float_t         refeta[99];   //[nref]
  Float_t         refy[99];   //[nref]
  Float_t         refphi[99];   //[nref]
  Float_t         refm[99];   //[nref]
  Float_t         refarea[99];   //[nref]
  Float_t         reftau1[99];   //[nref]
  Float_t         reftau2[99];   //[nref]
  Float_t         reftau3[99];   //[nref]
  Float_t         refdphijt[99];   //[nref]
  Float_t         refdrjt[99];   //[nref]
  Float_t         refparton_pt[99];   //[nref]
  Int_t           refparton_flavor[99];   //[nref]
  Int_t           refparton_flavorForB[99];   //[nref]
  Float_t         genChargedSum[99];   //[nref]
  Float_t         genHardSum[99];   //[nref]
  Float_t         signalChargedSum[99];   //[nref]
  Float_t         signalHardSum[99];   //[nref]
  Int_t           subid[99];   //[nref]
  Int_t           ngen;
  Int_t           genmatchindex[99];   //[ngen]
  Float_t         genpt[99];   //[ngen]
  Float_t         geneta[99];   //[ngen]
  Float_t         geny[99];   //[ngen]
  Float_t         gentau1[99];   //[ngen]
  Float_t         gentau2[99];   //[ngen]
  Float_t         gentau3[99];   //[ngen]
  Float_t         genphi[99];   //[ngen]
  Float_t         genm[99];   //[ngen]
  Float_t         gendphijt[99];   //[ngen]
  Float_t         gendrjt[99];   //[ngen]
  Int_t           gensubid[99];   //[ngen]

  // Declaration of leaf types: skim
  Int_t           pPAprimaryVertexFilter;
  Int_t           pBeamScrapingFilter;



  // List of branches: onia
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_QQ_type;   //!
  TBranch        *b_Reco_QQ_sign;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_mupl_idx;   //!
  TBranch        *b_Reco_QQ_mumi_idx;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_isCowboy;   //!
  TBranch        *b_Reco_QQ_ctau;   //!
  TBranch        *b_Reco_QQ_ctauErr;   //!
  TBranch        *b_Reco_QQ_cosAlpha;   //!
  TBranch        *b_Reco_QQ_ctau3D;   //!
  TBranch        *b_Reco_QQ_ctauErr3D;   //!
  TBranch        *b_Reco_QQ_cosAlpha3D;   //!
  TBranch        *b_Reco_QQ_whichGen;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!
  TBranch        *b_Reco_QQ_dca;   //!
  TBranch        *b_Reco_QQ_MassErr;   //!
  TBranch        *b_Reco_QQ_vtx;   //!
  TBranch        *b_Reco_QQ_Ntrk;   //!
  TBranch        *b_Reco_QQ_mupl_dxy_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mumi_dxy_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mupl_dxyErr_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mumi_dxyErr_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mupl_dz_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mumi_dz_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mupl_dzErr_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mumi_dzErr_muonlessVtx;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_mu_type;   //!
  TBranch        *b_Reco_mu_whichGen;   //!
  TBranch        *b_Reco_mu_SelectionType;   //!
  TBranch        *b_Reco_mu_charge;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_mu_trig;   //!
  TBranch        *b_Reco_mu_highPurity;   //!
  TBranch        *b_Reco_mu_TrkMuArb;   //!
  TBranch        *b_Reco_mu_TMOneStaTight;   //!
  TBranch        *b_Reco_mu_nPixValHits;   //!
  TBranch        *b_Reco_mu_nMuValHits;   //!
  TBranch        *b_Reco_mu_nTrkHits;   //!
  TBranch        *b_Reco_mu_normChi2_inner;   //!
  TBranch        *b_Reco_mu_normChi2_global;   //!
  TBranch        *b_Reco_mu_nPixWMea;   //!
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  TBranch        *b_Reco_mu_StationsMatched;   //!
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  TBranch        *b_Reco_mu_pt_inner;   //!
  TBranch        *b_Reco_mu_pt_global;   //!
  TBranch        *b_Reco_mu_ptErr_inner;   //!
  TBranch        *b_Reco_mu_ptErr_global;   //!
  TBranch        *b_Gen_weight;   //!
  TBranch        *b_Gen_pthat;   //!
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_QQ_type;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_QQ_momId;   //!
  TBranch        *b_Gen_QQ_ctau;   //!
  TBranch        *b_Gen_QQ_ctau3D;   //!
  TBranch        *b_Gen_QQ_mupl_idx;   //!
  TBranch        *b_Gen_QQ_mumi_idx;   //!
  TBranch        *b_Gen_QQ_whichRec;   //!
  TBranch        *b_Gen_mu_size;   //!
  TBranch        *b_Gen_mu_type;   //!
  TBranch        *b_Gen_mu_charge;   //!
  TBranch        *b_Gen_mu_4mom;   //!
  TBranch        *b_Gen_mu_whichRec;   //!

  // List of branches: jet
  TBranch        *b_evt;   //!
  TBranch        *b_b;   //!
  TBranch        *b_vx;   //!
  TBranch        *b_vy;   //!
  TBranch        *b_vz;   //!
  TBranch        *b_nref;   //!
  TBranch        *b_rawpt;   //!
  TBranch        *b_jtpt;   //!
  TBranch        *b_jteta;   //!
  TBranch        *b_jty;   //!
  TBranch        *b_jtphi;   //!
  TBranch        *b_jtpu;   //!
  TBranch        *b_jtm;   //!
  TBranch        *b_jtarea;   //!
  TBranch        *b_jtPfCHF;   //!
  TBranch        *b_jtPfNHF;   //!
  TBranch        *b_jtPfCEF;   //!
  TBranch        *b_jtPfNEF;   //!
  TBranch        *b_jtPfMUF;   //!
  TBranch        *b_jtPfCHM;   //!
  TBranch        *b_jtPfNHM;   //!
  TBranch        *b_jtPfCEM;   //!
  TBranch        *b_jtPfNEM;   //!
  TBranch        *b_jtPfMUM;   //!
  TBranch        *b_jttau1;   //!
  TBranch        *b_jttau2;   //!
  TBranch        *b_jttau3;   //!
  TBranch        *b_discr_jetID_cuts;   //!
  TBranch        *b_discr_jetID_bdt;   //!
  TBranch        *b_discr_fr01;   //!
  TBranch        *b_trackMax;   //!
  TBranch        *b_trackSum;   //!
  TBranch        *b_trackN;   //!
  TBranch        *b_trackHardSum;   //!
  TBranch        *b_trackHardN;   //!
  TBranch        *b_chargedMax;   //!
  TBranch        *b_chargedSum;   //!
  TBranch        *b_chargedN;   //!
  TBranch        *b_chargedHardSum;   //!
  TBranch        *b_chargedHardN;   //!
  TBranch        *b_photonMax;   //!
  TBranch        *b_photonSum;   //!
  TBranch        *b_photonN;   //!
  TBranch        *b_photonHardSum;   //!
  TBranch        *b_photonHardN;   //!
  TBranch        *b_neutralMax;   //!
  TBranch        *b_neutralSum;   //!
  TBranch        *b_neutralN;   //!
  TBranch        *b_hcalSum;   //!
  TBranch        *b_ecalSum;   //!
  TBranch        *b_eMax;   //!
  TBranch        *b_eSum;   //!
  TBranch        *b_eN;   //!
  TBranch        *b_muMax;   //!
  TBranch        *b_muSum;   //!
  TBranch        *b_muN;   //!
  TBranch        *b_beamId1;   //!
  TBranch        *b_beamId2;   //!
  TBranch        *b_pthat;   //!
  TBranch        *b_refpt;   //!
  TBranch        *b_refeta;   //!
  TBranch        *b_refy;   //!
  TBranch        *b_refphi;   //!
  TBranch        *b_refm;   //!
  TBranch        *b_refarea;   //!
  TBranch        *b_reftau1;   //!
  TBranch        *b_reftau2;   //!
  TBranch        *b_reftau3;   //!
  TBranch        *b_refdphijt;   //!
  TBranch        *b_refdrjt;   //!
  TBranch        *b_refparton_pt;   //!
  TBranch        *b_refparton_flavor;   //!
  TBranch        *b_refparton_flavorForB;   //!
  TBranch        *b_genChargedSum;   //!
  TBranch        *b_genHardSum;   //!
  TBranch        *b_signalChargedSum;   //!
  TBranch        *b_signalHardSum;   //!
  TBranch        *b_subid;   //!
  TBranch        *b_ngen;   //!
  TBranch        *b_genmatchindex;   //!
  TBranch        *b_genpt;   //!
  TBranch        *b_geneta;   //!
  TBranch        *b_geny;   //!
  TBranch        *b_gentau1;   //!
  TBranch        *b_gentau2;   //!
  TBranch        *b_gentau3;   //!
  TBranch        *b_genphi;   //!
  TBranch        *b_genm;   //!
  TBranch        *b_gendphijt;   //!
  TBranch        *b_gendrjt;   //!
  TBranch        *b_gensubid;   //!

  // List of branches: skim
  TBranch        *b_pPAprimaryVertexFilter;   //!
  TBranch        *b_pBeamScrapingFilter;   //!

  myTree(Bool_t ismc = false, Int_t nState=1, Bool_t jetR0p4 = true);
  virtual ~myTree();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     makeFlatTree();
  virtual Bool_t   areMuonsInAcceptance2018(Int_t iRecoQQ);
  virtual Bool_t   isGlobalMuonInAccept2018(TLorentzVector* Muon);
  virtual Bool_t   passQualityCuts2018(Int_t iRecoQQ);
  virtual Bool_t   isTriggerMatch(Int_t iRecoQQ, Int_t TriggerBit);
  virtual Bool_t   isMatchedDiMuon(Int_t iRecoQQ);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myTree_cxx
myTree::myTree(Bool_t ismc, Int_t nState, Bool_t jetR0p4) : fChain(0)  {
  isMC = ismc;
  jetR04 = jetR0p4;
  nS = nState;

  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  TTree *oniaTree = NULL;
  TTree *jetTree = NULL;
  TTree* skimTree = NULL;

  string fileName = "/data_CMS/cms/mnguyen/upsJet/data/merged_HiForestAOD.root";
  if (isMC) fileName = Form("/data_CMS/cms/mnguyen/upsJet/mc/%dS/merged_HiForestAOD.root",nS);


  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName.c_str());
  if (!f || !f->IsOpen()) {
    f = new TFile(fileName.c_str());
  }
  TDirectory * oniaDir = (TDirectory*)f->Get(Form("%s:/hionia",fileName.c_str()));
  oniaDir->GetObject("myTree",oniaTree);
  
  TDirectory * jetDir = (TDirectory*)f->Get(Form("%s:/ak%dPFJetAnalyzer",fileName.c_str(), jetR0p4?4:3));
  jetDir->GetObject("t",jetTree);

  TDirectory * skimDir = (TDirectory*)f->Get(Form("%s:/skimanalysis",fileName.c_str()));
  skimDir->GetObject("HltTree",skimTree);
  
  oniaTree->AddFriend(jetTree);
  oniaTree->AddFriend(skimTree);
  Init(oniaTree);
}

myTree::~myTree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t myTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t myTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void myTree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  Reco_QQ_4mom = 0;
  Reco_QQ_vtx = 0;
  Reco_mu_4mom = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  //onia branches
  if (fChain->GetBranch("Reco_QQ_size")) fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  if (fChain->GetBranch("Reco_QQ_type")) fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
  if (fChain->GetBranch("Reco_QQ_sign")) fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  if (fChain->GetBranch("Reco_QQ_4mom")) fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  if (fChain->GetBranch("Reco_QQ_mupl_idx")) fChain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
  if (fChain->GetBranch("Reco_QQ_mumi_idx")) fChain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
  if (fChain->GetBranch("Reco_QQ_trig")) fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  if (fChain->GetBranch("Reco_QQ_isCowboy")) fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
  if (fChain->GetBranch("Reco_QQ_ctau")) fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
  if (fChain->GetBranch("Reco_QQ_ctauErr")) fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
  if (fChain->GetBranch("Reco_QQ_cosAlpha")) fChain->SetBranchAddress("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha, &b_Reco_QQ_cosAlpha);
  if (fChain->GetBranch("Reco_QQ_ctau3D")) fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  if (fChain->GetBranch("Reco_QQ_ctauErr3D")) fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
  if (fChain->GetBranch("Reco_QQ_cosAlpha3D")) fChain->SetBranchAddress("Reco_QQ_cosAlpha3D", Reco_QQ_cosAlpha3D, &b_Reco_QQ_cosAlpha3D);
  if (isMC) {
    if (fChain->GetBranch("Reco_QQ_whichGen")) fChain->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen, &b_Reco_QQ_whichGen);
  }
  if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  if (fChain->GetBranch("Reco_QQ_dca")) fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
  if (fChain->GetBranch("Reco_QQ_MassErr")) fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
  if (fChain->GetBranch("Reco_QQ_vtx")) fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
  if (fChain->GetBranch("Reco_QQ_Ntrk")) fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
  if (fChain->GetBranch("Reco_QQ_mupl_dxy_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dxy_muonlessVtx", Reco_QQ_mupl_dxy_muonlessVtx, &b_Reco_QQ_mupl_dxy_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dxy_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dxy_muonlessVtx", Reco_QQ_mumi_dxy_muonlessVtx, &b_Reco_QQ_mumi_dxy_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mupl_dxyErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dxyErr_muonlessVtx", Reco_QQ_mupl_dxyErr_muonlessVtx, &b_Reco_QQ_mupl_dxyErr_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dxyErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dxyErr_muonlessVtx", Reco_QQ_mumi_dxyErr_muonlessVtx, &b_Reco_QQ_mumi_dxyErr_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mupl_dz_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dz_muonlessVtx", Reco_QQ_mupl_dz_muonlessVtx, &b_Reco_QQ_mupl_dz_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dz_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dz_muonlessVtx", Reco_QQ_mumi_dz_muonlessVtx, &b_Reco_QQ_mumi_dz_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mupl_dzErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dzErr_muonlessVtx", Reco_QQ_mupl_dzErr_muonlessVtx, &b_Reco_QQ_mupl_dzErr_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dzErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dzErr_muonlessVtx", Reco_QQ_mumi_dzErr_muonlessVtx, &b_Reco_QQ_mumi_dzErr_muonlessVtx);
  if (fChain->GetBranch("Reco_mu_size")) fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  if (fChain->GetBranch("Reco_mu_type")) fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
  if (isMC) {
    if (fChain->GetBranch("Reco_mu_whichGen")) fChain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
  }
  if (fChain->GetBranch("Reco_mu_SelectionType")) fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
  if (fChain->GetBranch("Reco_mu_charge")) fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
  if (fChain->GetBranch("Reco_mu_4mom")) fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  if (fChain->GetBranch("Reco_mu_trig")) fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  if (fChain->GetBranch("Reco_mu_highPurity")) fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
  if (fChain->GetBranch("Reco_mu_TrkMuArb")) fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
  if (fChain->GetBranch("Reco_mu_TMOneStaTight")) fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  if (fChain->GetBranch("Reco_mu_nPixValHits")) fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  if (fChain->GetBranch("Reco_mu_nMuValHits")) fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  if (fChain->GetBranch("Reco_mu_nTrkHits")) fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  if (fChain->GetBranch("Reco_mu_normChi2_inner")) fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
  if (fChain->GetBranch("Reco_mu_normChi2_global")) fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  if (fChain->GetBranch("Reco_mu_nPixWMea")) fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  if (fChain->GetBranch("Reco_mu_nTrkWMea")) fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  if (fChain->GetBranch("Reco_mu_StationsMatched")) fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  if (fChain->GetBranch("Reco_mu_dxy")) fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  if (fChain->GetBranch("Reco_mu_dxyErr")) fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  if (fChain->GetBranch("Reco_mu_dz")) fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  if (fChain->GetBranch("Reco_mu_dzErr")) fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  if (fChain->GetBranch("Reco_mu_pt_inner")) fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
  if (fChain->GetBranch("Reco_mu_pt_global")) fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
  if (fChain->GetBranch("Reco_mu_ptErr_inner")) fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
  if (fChain->GetBranch("Reco_mu_ptErr_global")) fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
  if (isMC) {
    if (fChain->GetBranch("Gen_weight")) fChain->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);
    if (fChain->GetBranch("Gen_pthat")) fChain->SetBranchAddress("Gen_pthat", &Gen_pthat, &b_Gen_pthat);
    if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
    if (fChain->GetBranch("Gen_QQ_type")) fChain->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
    if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
    if (fChain->GetBranch("Gen_QQ_momId")) fChain->SetBranchAddress("Gen_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
    if (fChain->GetBranch("Gen_QQ_ctau")) fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
    if (fChain->GetBranch("Gen_QQ_ctau3D")) fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
    if (fChain->GetBranch("Gen_QQ_mupl_idx")) fChain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
    if (fChain->GetBranch("Gen_QQ_mumi_idx")) fChain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
    if (fChain->GetBranch("Gen_QQ_whichRec")) fChain->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec, &b_Gen_QQ_whichRec);
    if (fChain->GetBranch("Gen_mu_size")) fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
    if (fChain->GetBranch("Gen_mu_type")) fChain->SetBranchAddress("Gen_mu_type", Gen_mu_type, &b_Gen_mu_type);
    if (fChain->GetBranch("Gen_mu_charge")) fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
    if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
    if (fChain->GetBranch("Gen_mu_whichRec")) fChain->SetBranchAddress("Gen_mu_whichRec", Gen_mu_whichRec, &b_Gen_mu_whichRec);
  }
   
   
  //jet branches
  if (fChain->GetBranch("evt")) fChain->SetBranchAddress("evt", &evt, &b_evt);
  if (fChain->GetBranch("b")) fChain->SetBranchAddress("b", &b, &b_b);
  if (fChain->GetBranch("vx")) fChain->SetBranchAddress("vx", &vx, &b_vx);
  if (fChain->GetBranch("vy")) fChain->SetBranchAddress("vy", &vy, &b_vy);
  if (fChain->GetBranch("vz")) fChain->SetBranchAddress("vz", &vz, &b_vz);
  if (fChain->GetBranch("nref")) fChain->SetBranchAddress("nref", &nref, &b_nref);
  if (fChain->GetBranch("rawpt")) fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
  if (fChain->GetBranch("jtpt")) fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
  if (fChain->GetBranch("jteta")) fChain->SetBranchAddress("jteta", jteta, &b_jteta);
  if (fChain->GetBranch("jty")) fChain->SetBranchAddress("jty", jty, &b_jty);
  if (fChain->GetBranch("jtphi")) fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
  if (fChain->GetBranch("jtpu")) fChain->SetBranchAddress("jtpu", jtpu, &b_jtpu);
  if (fChain->GetBranch("jtm")) fChain->SetBranchAddress("jtm", jtm, &b_jtm);
  if (fChain->GetBranch("jtarea")) fChain->SetBranchAddress("jtarea", jtarea, &b_jtarea);
  if (fChain->GetBranch("jtPfCHF")) fChain->SetBranchAddress("jtPfCHF", jtPfCHF, &b_jtPfCHF);
  if (fChain->GetBranch("jtPfNHF")) fChain->SetBranchAddress("jtPfNHF", jtPfNHF, &b_jtPfNHF);
  if (fChain->GetBranch("jtPfCEF")) fChain->SetBranchAddress("jtPfCEF", jtPfCEF, &b_jtPfCEF);
  if (fChain->GetBranch("jtPfNEF")) fChain->SetBranchAddress("jtPfNEF", jtPfNEF, &b_jtPfNEF);
  if (fChain->GetBranch("jtPfMUF")) fChain->SetBranchAddress("jtPfMUF", jtPfMUF, &b_jtPfMUF);
  if (fChain->GetBranch("jtPfCHM")) fChain->SetBranchAddress("jtPfCHM", jtPfCHM, &b_jtPfCHM);
  if (fChain->GetBranch("jtPfNHM")) fChain->SetBranchAddress("jtPfNHM", jtPfNHM, &b_jtPfNHM);
  if (fChain->GetBranch("jtPfCEM")) fChain->SetBranchAddress("jtPfCEM", jtPfCEM, &b_jtPfCEM);
  if (fChain->GetBranch("jtPfNEM")) fChain->SetBranchAddress("jtPfNEM", jtPfNEM, &b_jtPfNEM);
  if (fChain->GetBranch("jtPfMUM")) fChain->SetBranchAddress("jtPfMUM", jtPfMUM, &b_jtPfMUM);
  if (fChain->GetBranch("jttau1")) fChain->SetBranchAddress("jttau1", jttau1, &b_jttau1);
  if (fChain->GetBranch("jttau2")) fChain->SetBranchAddress("jttau2", jttau2, &b_jttau2);
  if (fChain->GetBranch("jttau3")) fChain->SetBranchAddress("jttau3", jttau3, &b_jttau3);
  if (fChain->GetBranch("discr_jetID_cuts")) fChain->SetBranchAddress("discr_jetID_cuts", discr_jetID_cuts, &b_discr_jetID_cuts);
  if (fChain->GetBranch("discr_jetID_bdt")) fChain->SetBranchAddress("discr_jetID_bdt", discr_jetID_bdt, &b_discr_jetID_bdt);
  if (fChain->GetBranch("discr_fr01")) fChain->SetBranchAddress("discr_fr01", discr_fr01, &b_discr_fr01);
  if (fChain->GetBranch("trackMax")) fChain->SetBranchAddress("trackMax", trackMax, &b_trackMax);
  if (fChain->GetBranch("trackSum")) fChain->SetBranchAddress("trackSum", trackSum, &b_trackSum);
  if (fChain->GetBranch("trackN")) fChain->SetBranchAddress("trackN", trackN, &b_trackN);
  if (fChain->GetBranch("trackHardSum")) fChain->SetBranchAddress("trackHardSum", trackHardSum, &b_trackHardSum);
  if (fChain->GetBranch("trackHardN")) fChain->SetBranchAddress("trackHardN", trackHardN, &b_trackHardN);
  if (fChain->GetBranch("chargedMax")) fChain->SetBranchAddress("chargedMax", chargedMax, &b_chargedMax);
  if (fChain->GetBranch("chargedSum")) fChain->SetBranchAddress("chargedSum", chargedSum, &b_chargedSum);
  if (fChain->GetBranch("chargedN")) fChain->SetBranchAddress("chargedN", chargedN, &b_chargedN);
  if (fChain->GetBranch("chargedHardSum")) fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
  if (fChain->GetBranch("chargedHardN")) fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
  if (fChain->GetBranch("photonMax")) fChain->SetBranchAddress("photonMax", photonMax, &b_photonMax);
  if (fChain->GetBranch("photonSum")) fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
  if (fChain->GetBranch("photonN")) fChain->SetBranchAddress("photonN", photonN, &b_photonN);
  if (fChain->GetBranch("photonHardSum")) fChain->SetBranchAddress("photonHardSum", photonHardSum, &b_photonHardSum);
  if (fChain->GetBranch("photonHardN")) fChain->SetBranchAddress("photonHardN", photonHardN, &b_photonHardN);
  if (fChain->GetBranch("neutralMax")) fChain->SetBranchAddress("neutralMax", neutralMax, &b_neutralMax);
  if (fChain->GetBranch("neutralSum")) fChain->SetBranchAddress("neutralSum", neutralSum, &b_neutralSum);
  if (fChain->GetBranch("neutralN")) fChain->SetBranchAddress("neutralN", neutralN, &b_neutralN);
  if (fChain->GetBranch("hcalSum")) fChain->SetBranchAddress("hcalSum", hcalSum, &b_hcalSum);
  if (fChain->GetBranch("ecalSum")) fChain->SetBranchAddress("ecalSum", ecalSum, &b_ecalSum);
  if (fChain->GetBranch("eMax")) fChain->SetBranchAddress("eMax", eMax, &b_eMax);
  if (fChain->GetBranch("eSum")) fChain->SetBranchAddress("eSum", eSum, &b_eSum);
  if (fChain->GetBranch("eN")) fChain->SetBranchAddress("eN", eN, &b_eN);
  if (fChain->GetBranch("muMax")) fChain->SetBranchAddress("muMax", muMax, &b_muMax);
  if (fChain->GetBranch("muSum")) fChain->SetBranchAddress("muSum", muSum, &b_muSum);
  if (fChain->GetBranch("muN")) fChain->SetBranchAddress("muN", muN, &b_muN);
  if (fChain->GetBranch("beamId1")) fChain->SetBranchAddress("beamId1", &beamId1, &b_beamId1);
  if (fChain->GetBranch("beamId2")) fChain->SetBranchAddress("beamId2", &beamId2, &b_beamId2);
  if (isMC) {
    if (fChain->GetBranch("pthat")) fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
    if (fChain->GetBranch("refpt")) fChain->SetBranchAddress("refpt", refpt, &b_refpt);
    if (fChain->GetBranch("refeta")) fChain->SetBranchAddress("refeta", refeta, &b_refeta);
    if (fChain->GetBranch("refy")) fChain->SetBranchAddress("refy", refy, &b_refy);
    if (fChain->GetBranch("refphi")) fChain->SetBranchAddress("refphi", refphi, &b_refphi);
    if (fChain->GetBranch("refm")) fChain->SetBranchAddress("refm", refm, &b_refm);
    if (fChain->GetBranch("refarea")) fChain->SetBranchAddress("refarea", refarea, &b_refarea);
    if (fChain->GetBranch("reftau1")) fChain->SetBranchAddress("reftau1", reftau1, &b_reftau1);
    if (fChain->GetBranch("reftau2")) fChain->SetBranchAddress("reftau2", reftau2, &b_reftau2);
    if (fChain->GetBranch("reftau3")) fChain->SetBranchAddress("reftau3", reftau3, &b_reftau3);
    if (fChain->GetBranch("refdphijt")) fChain->SetBranchAddress("refdphijt", refdphijt, &b_refdphijt);
    if (fChain->GetBranch("refdrjt")) fChain->SetBranchAddress("refdrjt", refdrjt, &b_refdrjt);
    if (fChain->GetBranch("refparton_pt")) fChain->SetBranchAddress("refparton_pt", refparton_pt, &b_refparton_pt);
    if (fChain->GetBranch("refparton_flavor")) fChain->SetBranchAddress("refparton_flavor", refparton_flavor, &b_refparton_flavor);
    if (fChain->GetBranch("refparton_flavorForB")) fChain->SetBranchAddress("refparton_flavorForB", refparton_flavorForB, &b_refparton_flavorForB);
    if (fChain->GetBranch("genChargedSum")) fChain->SetBranchAddress("genChargedSum", genChargedSum, &b_genChargedSum);
    if (fChain->GetBranch("genHardSum")) fChain->SetBranchAddress("genHardSum", genHardSum, &b_genHardSum);
    if (fChain->GetBranch("signalChargedSum")) fChain->SetBranchAddress("signalChargedSum", signalChargedSum, &b_signalChargedSum);
    if (fChain->GetBranch("signalHardSum")) fChain->SetBranchAddress("signalHardSum", signalHardSum, &b_signalHardSum);
    if (fChain->GetBranch("subid")) fChain->SetBranchAddress("subid", subid, &b_subid);
    if (fChain->GetBranch("ngen")) fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
    if (fChain->GetBranch("genmatchindex")) fChain->SetBranchAddress("genmatchindex", genmatchindex, &b_genmatchindex);
    if (fChain->GetBranch("genpt")) fChain->SetBranchAddress("genpt", genpt, &b_genpt);
    if (fChain->GetBranch("geneta")) fChain->SetBranchAddress("geneta", geneta, &b_geneta);
    if (fChain->GetBranch("geny")) fChain->SetBranchAddress("geny", geny, &b_geny);
    if (fChain->GetBranch("gentau1")) fChain->SetBranchAddress("gentau1", gentau1, &b_gentau1);
    if (fChain->GetBranch("gentau2")) fChain->SetBranchAddress("gentau2", gentau2, &b_gentau2);
    if (fChain->GetBranch("gentau3")) fChain->SetBranchAddress("gentau3", gentau3, &b_gentau3);
    if (fChain->GetBranch("genphi")) fChain->SetBranchAddress("genphi", genphi, &b_genphi);
    if (fChain->GetBranch("genm")) fChain->SetBranchAddress("genm", genm, &b_genm);
    if (fChain->GetBranch("gendphijt")) fChain->SetBranchAddress("gendphijt", gendphijt, &b_gendphijt);
    if (fChain->GetBranch("gendrjt")) fChain->SetBranchAddress("gendrjt", gendrjt, &b_gendrjt);
    if (fChain->GetBranch("gensubid")) fChain->SetBranchAddress("gensubid", gensubid, &b_gensubid);
  }

  if (fChain->GetBranch("pPAprimaryVertexFilter")) fChain->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter, &b_pPAprimaryVertexFilter);
  if (fChain->GetBranch("pBeamScrapingFilter")) fChain->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter, &b_pBeamScrapingFilter);

  fChain->SetBranchStatus("*",0);
  if (fChain->GetBranch("evt")) fChain->SetBranchStatus("evt",1);

  if (fChain->GetBranch("HLTriggers")) fChain->SetBranchStatus("HLTriggers",1);
  if (fChain->GetBranch("Reco_QQ_trig")) fChain->SetBranchStatus("Reco_QQ_trig",1);
  if (fChain->GetBranch("Reco_mu_trig")) fChain->SetBranchStatus("Reco_mu_trig",1);
  if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchStatus("Reco_QQ_VtxProb",1);
  if (fChain->GetBranch("Reco_mu_SelectionType")) fChain->SetBranchStatus("Reco_mu_SelectionType",1);
  if (fChain->GetBranch("Reco_mu_nTrkWMea")) fChain->SetBranchStatus("Reco_mu_nTrkWMea",1);
  if (fChain->GetBranch("Reco_mu_nPixWMea")) fChain->SetBranchStatus("Reco_mu_nPixWMea",1);
  if (fChain->GetBranch("Reco_mu_dxy")) fChain->SetBranchStatus("Reco_mu_dxy",1);
  if (fChain->GetBranch("Reco_mu_dz")) fChain->SetBranchStatus("Reco_mu_dz",1);
  if (fChain->GetBranch("Reco_mu_nTrkHits")) fChain->SetBranchStatus("Reco_mu_nTrkHits",1);
  if (fChain->GetBranch("Reco_mu_normChi2_global")) fChain->SetBranchStatus("Reco_mu_normChi2_global",1);
  if (fChain->GetBranch("Reco_mu_normChi2_inner")) fChain->SetBranchStatus("Reco_mu_normChi2_inner",1);
  if (fChain->GetBranch("Reco_mu_TrkMuArb")) fChain->SetBranchStatus("Reco_mu_TrkMuArb",1);
  if (fChain->GetBranch("Reco_QQ_size")) fChain->SetBranchStatus("Reco_QQ_size",1);
  if (fChain->GetBranch("Reco_QQ_sign")) fChain->SetBranchStatus("Reco_QQ_sign",1);
  if (fChain->GetBranch("Reco_QQ_4mom")) fChain->SetBranchStatus("Reco_QQ_4mom",1);
  if (fChain->GetBranch("Reco_QQ_mupl_idx")) fChain->SetBranchStatus("Reco_QQ_mupl_idx",1);
  if (fChain->GetBranch("Reco_QQ_mumi_idx")) fChain->SetBranchStatus("Reco_QQ_mumi_idx",1);
  if (fChain->GetBranch("Reco_mu_4mom")) fChain->SetBranchStatus("Reco_mu_4mom",1);
  if (fChain->GetBranch("Reco_mu_size")) fChain->SetBranchStatus("Reco_mu_size",1);
  if (isMC) {
    if (fChain->GetBranch("Reco_QQ_whichGen")) fChain->SetBranchStatus("Reco_QQ_whichGen",1);
    if (fChain->GetBranch("Reco_mu_whichGen")) fChain->SetBranchStatus("Reco_mu_whichGen",1);
    if (fChain->GetBranch("Gen_weight")) fChain->SetBranchStatus("Gen_weight",1);
    if (fChain->GetBranch("Gen_pthat")) fChain->SetBranchStatus("Gen_pthat",1);
    if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchStatus("Gen_QQ_4mom",1);
    if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchStatus("Gen_QQ_size",1);
    if (fChain->GetBranch("Gen_QQ_mupl_idx")) fChain->SetBranchStatus("Gen_QQ_mupl_idx",1);
    if (fChain->GetBranch("Gen_QQ_mumi_idx")) fChain->SetBranchStatus("Gen_QQ_mumi_idx",1);
    if (fChain->GetBranch("Gen_QQ_ctau3D")) fChain->SetBranchStatus("Gen_QQ_ctau3D",1);
    if (fChain->GetBranch("Gen_QQ_ctau")) fChain->SetBranchStatus("Gen_QQ_ctau",1);
    if (fChain->GetBranch("Gen_QQ_whichRec")) fChain->SetBranchStatus("Gen_QQ_whichRec",1);
    if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchStatus("Gen_mu_4mom",1);
    if (fChain->GetBranch("Gen_mu_size")) fChain->SetBranchStatus("Gen_mu_size",1);
    if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchStatus("Gen_QQ_mupl_4mom",1);
    if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchStatus("Gen_QQ_mumi_4mom",1);
  }

  if (fChain->GetBranch("nref")) fChain->SetBranchStatus("nref",1);
  if (fChain->GetBranch("rawpt")) fChain->SetBranchStatus("rawpt",1);
  if (fChain->GetBranch("jtpt")) fChain->SetBranchStatus("jtpt",1);
  if (fChain->GetBranch("jteta")) fChain->SetBranchStatus("jteta",1);
  if (fChain->GetBranch("jty")) fChain->SetBranchStatus("jty",1);
  if (fChain->GetBranch("jtphi")) fChain->SetBranchStatus("jtphi",1);
  if (fChain->GetBranch("jtm")) fChain->SetBranchStatus("jtm",1);
  if (isMC) {
    if (fChain->GetBranch("pthat")) fChain->SetBranchStatus("pthat",1);
    if (fChain->GetBranch("refpt")) fChain->SetBranchStatus("refpt",1);
    if (fChain->GetBranch("refeta")) fChain->SetBranchStatus("refeta",1);
    if (fChain->GetBranch("refy")) fChain->SetBranchStatus("refy",1);
    if (fChain->GetBranch("refphi")) fChain->SetBranchStatus("refphi",1);
    if (fChain->GetBranch("refm")) fChain->SetBranchStatus("refm",1);
    if (fChain->GetBranch("refarea")) fChain->SetBranchStatus("refarea",1);
  }

  if (fChain->GetBranch("pPAprimaryVertexFilter")) fChain->SetBranchStatus("pPAprimaryVertexFilter",1);
  if (fChain->GetBranch("pBeamScrapingFilter")) fChain->SetBranchStatus("pBeamScrapingFilter",1);

  Notify();
}

Bool_t myTree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void myTree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t myTree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

Bool_t myTree::areMuonsInAcceptance2018(Int_t iRecoQQ){
  TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[iRecoQQ]);
  TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[iRecoQQ]);

  return ( isGlobalMuonInAccept2018(RecoQQmupl) && isGlobalMuonInAccept2018(RecoQQmumi) );
}

Bool_t myTree::isGlobalMuonInAccept2018(TLorentzVector* Muon) {
  //return (fabs(Muon->Eta()) < 2.4 &&
  //	  ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
  //	   (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.47-1.89*fabs(Muon->Eta())) ||
  //	   (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.5)));
  if (Muon->Pt() < 3.5) return false;
  return true;

}

Bool_t myTree::passQualityCuts2018(Int_t iRecoQQ){
  int iMupl = Reco_QQ_mupl_idx[iRecoQQ];
  int iMumi = Reco_QQ_mumi_idx[iRecoQQ];

  if ( ! (Reco_mu_SelectionType[iMumi]&((ULong64_t)pow(2, 1))) ) return false; // require the muons to be global muons
  if ( ! (Reco_mu_SelectionType[iMumi]&((ULong64_t)pow(2, 3))) ) return false; // require the muons to be tracker muons
  // if ( ! (Reco_mu_highPurity[iMumi]) ) return false;
  //if ( ! (Reco_mu_TMOneStaTight[iMumi]==1) ) return false; // = isGoodMuon
  if ( ! (Reco_mu_nTrkWMea[iMumi] > 5) ) return false;
  if ( ! (Reco_mu_nPixWMea[iMumi] > 0) ) return false;
  if ( ! (fabs(Reco_mu_dxy[iMumi]) < 0.3) ) return false;
  if ( ! (fabs(Reco_mu_dz[iMumi]) < 20.) ) return false;
    
  if ( ! (Reco_mu_SelectionType[iMupl]&((ULong64_t)pow(2, 1))) ) return false; // require the muons to be global muons
  if ( ! (Reco_mu_SelectionType[iMupl]&((ULong64_t)pow(2, 3))) ) return false; // require the muoons to be tracker muons
  // if ( ! (Reco_mu_highPurity[iMupl]) ) return false;
  //if ( ! (Reco_mu_TMOneStaTight[iMupl]==1) ) return false; // = isGoodMuon
  if ( ! (Reco_mu_nTrkWMea[iMupl] > 5) ) return false;
  if ( ! (Reco_mu_nPixWMea[iMupl] > 0) ) return false;
  if ( ! (fabs(Reco_mu_dxy[iMupl]) < 0.3) ) return false;
  if ( ! (fabs(Reco_mu_dz[iMupl]) < 20.) ) return false;
    
  if ( ! (Reco_QQ_VtxProb[iRecoQQ] > 0.01) ) return false;

  return true;
}

Bool_t myTree::isTriggerMatch(Int_t iRecoQQ, Int_t TriggerBit){
  
  //if (!( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) )) return false;
  if (!( (Reco_QQ_trig[iRecoQQ]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) )) return false;

  return true;
}

Bool_t myTree::isMatchedDiMuon(Int_t iRecoQQ){
  if (Reco_QQ_whichGen[iRecoQQ]<0) return false;
  if (Reco_mu_whichGen[Reco_QQ_mupl_idx[iRecoQQ]]<0) return false;
  if (Reco_mu_whichGen[Reco_QQ_mumi_idx[iRecoQQ]]<0) return false;
  return true;
}

double jecCorr(double jtPt, double rawPt, double jpsiPt) {
  return ( (1-(jpsiPt/rawPt))*jtPt + ((jpsiPt/rawPt)/(jtPt/rawPt))*jtPt );
  //return ( (1-z)*jtPt + z*rawPt );
}

double zjecCorr(double jtPt, double rawPt, double z) { //same as the one before but expressed in another way for a test
  //return ( (1-(jpsiPt/rawPt))*jtPt + ((jpsiPt/rawPt)/(jtPt/rawPt))*jtPt );
  return ( (1-z)*jtPt + z*rawPt );
}

double jeuCorr (double jtPt, double z, double jeu) {
  return ( (1-z)*(1+jeu)*jtPt + z*jtPt );
}

double getAccEff (double rap, double pt) {
  // reading the efficency from histograms
  return 1;
}

#endif // #ifdef myTree_cxx
