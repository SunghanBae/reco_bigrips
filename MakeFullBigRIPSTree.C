
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtRecoPID.hh"
#include "TArtCalibPID.hh"
#include "TArtCalibCoin.hh"
//#include "TArtCalibTSRef.hh"
#include "TArtEventInfo.hh"


#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"

#include "signal.h"

#include <string>
#include <iostream>

using namespace std;

// function to exit loop at keyboard interrupt. 
bool stoploop = false;
void stop_interrupt()
{
  printf("keyboard interrupt\n");
  stoploop = true;
}

//void MakeFullBigRIPSTree(char *infile, char *outfile="bigrips.root"){
int main(int argc, char** argv){

  //  signal(SIGINT,stop_interrupt); // CTRL + C , interrupt

  //  gSystem->Load("libXMLParser.so");
  //  gSystem->Load("libanaroot.so");

  TArtStoreManager * sman = TArtStoreManager::Instance();
  TArtEventStore *estore = new TArtEventStore();
  estore->SetInterrupt(&stoploop);
//  estore->Open(infile);
  estore->Open(argv[1]);
  //  estore->Open(2); // bigrips online

  TArtBigRIPSParameters *para = TArtBigRIPSParameters::Instance();
  para->LoadParameter("Fromkathrin_db/BigRIPSPPAC_1465435445_2016.xml");
  para->LoadParameter("db/ppacdefault.xml");
  para->LoadParameter("Fromkathrin_db/BigRIPSPlastic_2016.xml");
  para->LoadParameter("Fromkathrin_db/BigRIPSIC_2016.xml");
  para->LoadParameter("Fromkathrin_db/FocalPlane_2016.xml");
//  para->LoadParameter("Fromkathrin_db/BigRIPSGe.xml");

  TArtCalibPID *brcalib = new TArtCalibPID();
//  TArtCalibCoin *calcoin= new TArtCalibCoin();
  TArtRecoPID  *brreco = new TArtRecoPID();

  // setting up rips parameters
//  TArtRIPS * rips3to5 = brreco->DefineNewRIPS(3,5,"matrix/mat1.mat",7.073); // f3 - f5
//  TArtRIPS * rips5to7 = brreco->DefineNewRIPS(5,7,"matrix/mat2.mat",6.2675); // f5 - f7

///  TArtRIPS * rips3to5 = brreco->DefineNewRIPS(3,5,"matrix/mat1_sumikama.mat",5.3833); // f3 - f5
  // TArtRIPS * rips5to7 = brreco->DefineNewRIPS(5,7,"matrix/mat2.mat",4.9469); // f5 - f7
  TArtRIPS * rips3to5 = brreco->DefineNewRIPS(3,5,"mymat/mat35.mat",7.7565); // f3 - f5
//  TArtRIPS * rips5to7 = brreco->DefineNewRIPS(5,7,"matrix/mat2.mat",7.749); // f5 - f7
  TArtRIPS * rips5to11 = brreco->DefineNewRIPS(5,11,"mymat/mat511.mat",7.749); // f5 - f11


  // setting up tof
  //TArtTOF * tof3to7  = brreco->DefineNewTOF("F3pl","F7pl",248.0,5);//171Dy setting
  TArtTOF * tof3to11  = brreco->DefineNewTOF("F3pl","F11pl-1",446.98,5);//2016_June_Estrade_Bae
  //TArtTOF * tof3to7_cfd  = brreco->DefineNewTOF("F3pl_cfd","F7pl_cfd",250.5,5);
  //  TArtTOF * tof3to7  = brreco->DefineNewTOF("F3pl","F7pl",275.733,5);

  // setting up beam id devices
//  TArtBeam *beam37ic7 = brreco->DefineNewBeam(rips3to5,rips5to7,tof3to7,"F7IC");
  //TArtBeam *beam37ic7_cfd = brreco->DefineNewBeam(rips3to5,rips5to7,tof3to7_cfd,"F7IC");
  TArtBeam *beam37ic11 = brreco->DefineNewBeam(rips3to5,rips5to11,tof3to11,"F11IC");

  //  TArtBeam *beam37ic11_cfd = brreco->DefineNewBeam(rips3to5,rips5to7,tof3to7_cfd,"F11IC");

//  TFile *fout = new TFile(outfile,"RECREATE");
  TFile *fout = new TFile(argv[2],"RECREATE");
  TTree *tree = new TTree("tree","tree");

  Int_t maxbit;
  tree->Branch("maxbit",&maxbit,"maxbit/I");
  // define data nodes which are supposed to be dumped to tree 
  TClonesArray * info_array = 
    (TClonesArray *)sman->FindDataContainer("EventInfo");
  tree->Branch(info_array->GetName(),&info_array);
  TClonesArray * ppac_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSPPAC");
  tree->Branch(ppac_array->GetName(),&ppac_array);
  TClonesArray * pla_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSPlastic");
  tree->Branch(pla_array->GetName(),&pla_array);
  TClonesArray * ic_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSIC");
  tree->Branch(ic_array->GetName(),&ic_array);
  TClonesArray * ge_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSGe");
  tree->Branch(ge_array->GetName(),&ge_array);
  TClonesArray * fpl_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSFocalPlane");
  tree->Branch(fpl_array->GetName(),&fpl_array);
  TClonesArray * tof_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSTOF");
  tree->Branch(tof_array->GetName(),&tof_array);
  TClonesArray * rips_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSRIPS");
  tree->Branch(rips_array->GetName(),&rips_array);
  TClonesArray * beam_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSBeam");
  tree->Branch(beam_array->GetName(),&beam_array);

  //  TArtBeam *beam = (TArtBeam*)beam_array->At(0);

  cout << "Before loop " << endl;
  int neve = 0;
  while(estore->GetNextEvent()){
    //  while(estore->GetNextEvent()&&neve<100000){
    maxbit = 0;
    if(neve%1000==0){
      cout << "\revent: " << neve;
      cout.flush();}
//    if (neve>4000) break;
//    calcoin->ClearData();
//    calcoin->LoadData();
    brcalib->ClearData();
    brcalib->ReconstructData();
    brreco->ClearData();
    brreco->ReconstructData();

    //    beam->SetAoQ(beam->GetAoQ()-0.042);

    Int_t tbit = ((TArtEventInfo*)info_array->At(0))->GetTriggerBit();
    for(int i=1;;i++){
      if(tbit%2 == 1)
	maxbit = i;
      else break;

      tbit = tbit/2;
    }

    tree->Fill();

    neve ++;
    //if (neve>1000000) break;
  }

  cout << "\revent: " << neve << endl;
  fout->Write();
  fout->Close();

  return 0;
}

