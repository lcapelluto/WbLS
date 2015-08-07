/* Plot data from offset investigations. Data from samples
 * are in box measurements. */
void Offset() {
    
    TCanvas *c3 = new TCanvas("c3","c3",600, 400);
    
    double x_popop[] = {384, 407, 430};
    double y_popop[] = {16, 16, 16};
    
    double x_anthra[] = {375,397,420,446}; 
    double y_anthra[] = {7, 5, 5, 5};

    double x_ptp[] = {323,338,352};
    double new_ref[] = {325,339};
    double new_ref_offset[] = {4,3};
    double y_ptp[] = {6, 4, 6};
    double x_ptp_err[] = {0,0,3};
    double y_ptp_err[] = {0,0,2};
    
    double x_ppo[] = {355,371};
    double y_ppo[] = {4,4};

    //double x_emissionb[] = {300,350,400,450,475,477,500,525,550,600,650};
    //double y_emissionb[] = {1, 1,1,1.5,2,2,1.5,1,1.5,1.25,1};
    
    double x_emissionb[] = {350,380,400,450,500,550,600};
    double y_emissionb[] = {0,-.5,0,-.5,0,-.5,0};

    //excitation scans where x is the emission wavelength collected
    //double x_emissionIS[] = {350,400,475,500,545,600,645};
    //double y_emissionIS[] = {2.5,2.5,3.5,2.5,3,2.5,0.5};
    double x_emissionIS[] = {350,380,400,450,500,550,600};
    double y_emissionIS[] = {-2,-2.5,-2,-2.5,-2,-2.5,-2};
    
    double x_excitationb[] = {300,350,400,450,475,500,550,600,650};
    double y_excitationb[] = {0,-1,-1,-1,-1.5,-1,-1.5,-.5,-1.5};
    
    double x_excitationIS[] = {345,400,450,475,500,550,600,650};
    double y_excitationIS[] = {-2,-2,-2.5,-3.5,-2.5,-2.5,-2.5,-3};
    
    double x_cwf[] = {365,404.5,436,488,545,587,599,612,631};
    double y_cwf[] = {-4,-3,-4,-3,-4,-3.5,-3,-3,-3.5};
//    double y_cwf[] = {-.5,1,.5,1.5,.5,1,1.5,1.5,.5};
    double err[] = {.5,.5,.5,.5,.5,.5,.5,.5,.5};
    
    double x_cwf_b[] = {404.6,436,488,545,587,599,612,631,650.5,663};
    double y_cwf_b[] = {-1.5,-2,-1.5,-2,-1,-1,-1,-1,-1,-1};
//    double y_cwf_b[] = {2.5,2,2,2.5,3,3,3,3,3,3,2};
    //584

    TGraph *g1 = new TGraph(3, x_popop, y_popop);
    g1->SetTitle("POPOP");
    g1->SetMarkerStyle(23);
    g1->SetMarkerColor(4);
    g1->SetFillColor(0);
    g1->SetMarkerSize(1.2);
    
    TGraph *g2 = new TGraph(4, x_anthra, y_anthra);
    g2->SetTitle("Anthracene");
    g2->SetMarkerStyle(22);
    g2->SetMarkerColor(9);
    g2->SetFillColor(0);
    g2->SetMarkerSize(1.2);

    TGraphErrors *g3 = new TGraphErrors(3, x_ptp, y_ptp,x_ptp_err,y_ptp_err);
    g3->SetTitle("P-TP");
    g3->SetMarkerStyle(23);
    g3->SetMarkerColor(7);
    g3->SetFillColor(0);
    g3->SetMarkerSize(1.2);
    
    TGraphErrors *g4 = new TGraphErrors(7, x_emissionb, y_emissionb);
    g4->SetTitle("Excitation Scan");
    g4->SetMarkerStyle(25);
    g4->SetMarkerColor(2);
    g4->SetFillColor(0);
    g4->SetMarkerSize(1.2);
    
    TGraphErrors *g5 = new TGraphErrors(7, x_emissionIS, y_emissionIS);
    g5->SetTitle("Excitation Scan in IS");
    g5->SetMarkerStyle(20);
    g5->SetMarkerColor(6);
    g5->SetFillColor(0);
    g5->SetMarkerSize(1.2);
    
    TGraph *g6 = new TGraph(9, x_excitationb, y_excitationb);
    g6->SetTitle("Emission Scan");
    g6->SetMarkerStyle(25);
    g6->SetMarkerColor(3);
    g6->SetFillColor(0);
    g6->SetMarkerSize(1.2);
    g6->GetYaxis()->SetRangeUser(-4,17);
    g6->GetYaxis()->SetTitle("offset from expected [nm]");
    g6->GetXaxis()->SetTitle("expected wavelength [nm]");
    
    TGraph *g7 = new TGraph(8, x_excitationIS, y_excitationIS);
    g7->SetTitle("Emission Scan in IS");
    g7->SetMarkerStyle(20);
    g7->SetMarkerColor(8);
    g7->SetFillColor(0);
    g7->SetMarkerSize(1.2);
      
    TGraph *g8 = new TGraph(2, x_ppo, y_ppo);
    g8->SetTitle("PPO");
    g8->SetMarkerStyle(37);
    g8->SetMarkerColor(kViolet+2);
    g8->SetFillColor(0);
    g8->SetMarkerSize(1.2);
    
    TGraphErrors *g9 = new TGraphErrors(9, x_cwf, y_cwf);
    g9->SetTitle("Fluorescent lights in IS");
    g9->SetMarkerStyle(20);
    g9->SetMarkerColor(1);
    g9->SetFillColor(0);
    g9->SetMarkerSize(1.2);
    
    TGraphErrors *g10 = new TGraphErrors(10, x_cwf_b, y_cwf_b);
    g10->SetTitle("Fluorescent lights");
    g10->SetMarkerStyle(25);
    g10->SetMarkerColor(1);
    g10->SetFillColor(0);
    g10->SetMarkerSize(1.2);
    
    g6->Draw();
    //g4->Draw("SAME");
    //g5->Draw("SAME");
    g7->Draw("SAME");
    //g1->Draw("SAME");
    //g2->Draw("SAME");
    //g3->Draw("SAME");
    //g8->Draw("SAME");
    g9->Draw("SAME");
    g10->Draw("SAME");
    c3->BuildLegend();
}
