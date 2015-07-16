/* Plot concentration against QY.
 */
void GraphConcn() {
    double x[] = {1.12,10.0,48.5,100.4,890.8};
    double y[] = {.181,.121,.143,.1309,.0831};
    double yerror[] = {.005,.001,.0006,.0006,.0004};

    TGraphErrors* plot = new TGraphErrors(5, x, y,NULL,yerror);

    plot->SetTitle("Anthracene in cyclohexane");
    plot->GetXaxis()->SetTitle("Concentration [mg/L]");
    plot->GetYaxis()->SetTitle("QY");
    plot->Draw();
}

