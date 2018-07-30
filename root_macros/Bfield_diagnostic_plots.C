void Bfield_diagnostic_plots(const char *rootfilename){
  TFile *fin = new TFile(rootfilename,"READ");

  TH2F *field_x,*field_y,*field_z, *field;

  fin->GetObject("field_x",field_x);
  fin->GetObject("field_y",field_y);
  fin->GetObject("field_z",field_z);
  fin->GetObject("field",field);

  TH2F *hBx_yzproj_Earm,*hBy_yzproj_Earm,*hBz_yzproj_Earm,*hBtot_yzproj_Earm;
  TH2F *hBx_yzproj_Harm,*hBy_yzproj_Harm,*hBz_yzproj_Harm,*hBtot_yzproj_Harm;

  fin->GetObject("hBx_yzproj_Earm", hBx_yzproj_Earm);
  fin->GetObject("hBy_yzproj_Earm", hBy_yzproj_Earm);
  fin->GetObject("hBz_yzproj_Earm", hBz_yzproj_Earm);
  fin->GetObject("hBtot_yzproj_Earm", hBtot_yzproj_Earm);

  fin->GetObject("hBx_yzproj_Harm", hBx_yzproj_Harm);
  fin->GetObject("hBy_yzproj_Harm", hBy_yzproj_Harm);
  fin->GetObject("hBz_yzproj_Harm", hBz_yzproj_Harm);
  fin->GetObject("hBtot_yzproj_Harm", hBtot_yzproj_Harm);
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->SetGrid();
  c1->Divide(2,2,.001,.001);

  c1->cd(1)->SetGrid();

  field_x->Draw("colz");

  c1->cd(2)->SetGrid();
  field_y->Draw("colz");

  c1->cd(3)->SetGrid();
  field_z->Draw("colz");

  c1->cd(4)->SetGrid();
  field->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
  c2->SetGrid();
  c2->Divide(2,2,.001,.001);

  c2->cd(1)->SetGrid();

  hBx_yzproj_Earm->Draw("colz");

  c2->cd(2)->SetGrid();
  hBy_yzproj_Earm->Draw("colz");

  c2->cd(3)->SetGrid();
  hBz_yzproj_Earm->Draw("colz");

  c2->cd(4)->SetGrid();
  hBtot_yzproj_Earm->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","c3",1200,1200);
  c3->SetGrid();
  c3->Divide(2,2,.001,.001);
  
  c3->cd(1)->SetGrid();

  hBx_yzproj_Harm->Draw("colz");

  c3->cd(2)->SetGrid();
  hBy_yzproj_Harm->Draw("colz");

  c3->cd(3)->SetGrid();
  hBz_yzproj_Harm->Draw("colz");

  c3->cd(4)->SetGrid();
  hBtot_yzproj_Harm->Draw("colz");

}
