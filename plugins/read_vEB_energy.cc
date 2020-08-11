#include "ProdTutorial/ProducerTest/plugins/DetImgProducer.h"
#include <iostream>
using namespace std;

std::vector<float>& DetImgProducer::read_vEB_energy(int vec_size)
 {
    std::vector<float> *vEB_energy_read = 0;
    TFile *fr = TFile::Open("ECAL_Rechit.root","READ");
    if (!fr) { return *vEB_energy_read; }
 
    
    TTree *tr = (TTree*)fr->Get("vec_tree");
    tr->SetBranchAddress("vEB_energy_vec",&vEB_energy_read);
    for (int i=0;i<vec_size;i++){
      tr->GetEntry(i);
    } 
    for (int j=0;j<10;j++){
     std::cout<<vEB_energy_read->at(j)<<" ";
    }
    std::cout<<endl;
    tr->ResetBranchAddresses();
    fr->Close();
    return *vEB_energy_read;
 }
