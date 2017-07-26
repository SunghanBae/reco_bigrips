
void toff(int run){

TString fname = Form("run%04d.root",run);

TFile* f = new TFile(fname);
TTree* t = (TTree*)f->Get("tree");
TH1F* h = new TH1F("h","",100,0,30);
TCut cc[10];
Double_t aoq[10],y,c,brhou,brhod,lu,ld,mu,e,betau,betad;
TGraphErrors* g = new TGraphErrors(3);
//run1179 setting and length(mm) from Kathrin's Focalplane.xml 
lu=54917-31633;
ld=125984-54917;
brhou=7.7565;
brhod=7.751533;
c=299792458;
e=1.602176565*1e-19;
mu=931.4940954*1e+6/c/c;

//cuts for reference run 1179, with Kathrin's xmls
/*
cc[0] = "BigRIPSIC.fRawADCSqSum[2]> 7900 && BigRIPSIC.fRawADCSqSum[2] <8400 && BigRIPSTOF.tof >442.5 &&BigRIPSTOF.tof<444.5"; //88Se
cc[1] = "BigRIPSIC.fRawADCSqSum[2]> 8000 && BigRIPSIC.fRawADCSqSum[2] <8500 && BigRIPSTOF.tof >444.8 &&BigRIPSTOF.tof<446.8"; //89Se
cc[2] = "BigRIPSIC.fRawADCSqSum[2]> 8100 && BigRIPSIC.fRawADCSqSum[2] <8500 && BigRIPSTOF.tof >447.6 &&BigRIPSTOF.tof<449.6"; //90Se
*/
/*
//cuts for reference run 1190, with Kathrin's xmls
cc[0] = "BigRIPSIC.fRawADCSqSum[2]> 8300 && BigRIPSIC.fRawADCSqSum[2] <8700 && BigRIPSTOF.tof >444.5 &&BigRIPSTOF.tof<446.5"; //88Se
cc[1] = "BigRIPSIC.fRawADCSqSum[2]> 8400 && BigRIPSIC.fRawADCSqSum[2] <8800 && BigRIPSTOF.tof >447.5 &&BigRIPSTOF.tof<449.5"; //89Se
cc[2] = "BigRIPSIC.fRawADCSqSum[2]> 8500 && BigRIPSIC.fRawADCSqSum[2] <8900 && BigRIPSTOF.tof >450 &&BigRIPSTOF.tof<452"; //90Se
*/
/*
//cuts for n-rich run 1181, with Kathrin's xmls
cc[0] = "BigRIPSIC.fRawADCSqSum[2]> 8150 && BigRIPSIC.fRawADCSqSum[2] <8500 && BigRIPSTOF.tof >445.5 &&BigRIPSTOF.tof<447"; //89Se
cc[1] = "BigRIPSIC.fRawADCSqSum[2]> 8250 && BigRIPSIC.fRawADCSqSum[2] <8600 && BigRIPSTOF.tof >447 &&BigRIPSTOF.tof<450"; //90Se
cc[2] = "BigRIPSIC.fRawADCSqSum[2]> 8350 && BigRIPSIC.fRawADCSqSum[2] <8700 && BigRIPSTOF.tof >450 &&BigRIPSTOF.tof<452.5"; //91Se
*/

//cuts for n-rich run 1193, with Kathrin's xmls
cc[0] = "BigRIPSIC.fRawADCSqSum[2]> 8450 && BigRIPSIC.fRawADCSqSum[2] <8800 && BigRIPSTOF.tof >448 && BigRIPSTOF.tof<449"; //89Se
cc[1] = "BigRIPSIC.fRawADCSqSum[2]> 8550 && BigRIPSIC.fRawADCSqSum[2] <8900 && BigRIPSTOF.tof >450 && BigRIPSTOF.tof<452"; //90Se
cc[2] = "BigRIPSIC.fRawADCSqSum[2]> 8650 && BigRIPSIC.fRawADCSqSum[2] <9000 && BigRIPSTOF.tof >452.8 && BigRIPSTOF.tof<454.5"; //91Se

//reference run
/*
aoq[0] = 88./34.; // 88Se
aoq[1] = 89./34.; // 89Se
aoq[2] = 90./34.; // 90Se
*/

//n-rich run

aoq[0] = 89./34.; // 89Se
aoq[1] = 90./34.; // 90Se
aoq[2] = 91./34.; // 91Se



for( int i =0; i<3; i++){
t->Draw("BigRIPSPlastic.fTime[9]-BigRIPSPlastic.fTime[1]>>h",cc[i]);
betau=1./ sqrt(1+pow(mu*c*aoq[i]/brhou,2));
betad=1./ sqrt(1+pow(mu*c*aoq[i]/brhod,2));
cout<<"beta35 : "<< betau<<" beta511 : "<<betad<<endl;
y=(lu/betau+ld/betad)/c*1000000;
g->SetPoint(i,h->GetMean(),y);
g->SetPointError(i,h->GetStdDev(),0.00001*y);
}
TFile*nf = new TFile(Form("run%04d_toff.root",run),"RECREATE");
g->Draw("AP");
nf->WriteTObject(g);
nf->Close();

}


