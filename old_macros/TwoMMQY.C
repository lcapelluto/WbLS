#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Array size for x and y points. */
const int MAX_NUM_POINTS = 500;
/* Estimated maximum length of data file line. */
const int LINE_SIZE = 60;
/* 'N'umber of 'A'veraging 'P'oints to establish a baseline. */
const int NAP = 20;

/* Apply the selected emission correction, subtract offset, graph
 * spectrum and print the 2 measurement method Quantum Yield calculation.
 */
void TwoMMQY(char* scan = "anthracene_ethanol/AE_340_IN_default",
           double L_a = 0.0,
           char* title = "Graph",
           char* corr_table = "defaultCorrectionData",
           int corrheader = 0, int inc = 10, int offset = 300,
           int lbmin = 333, int lbmax = 344, int pbmin = 376, int pbmax = 438,
           int num_points = 351, int dataheader = 4,
           char* xtitle = "wavelength [nm]",
           char* ytitle = "Intensity [photons/sec]") {

    double QY;
    double L_b = 0.0, P_b = 0.0;
    double baseline = 0.0;
    double laerror = 33868.729016;
    double lberror = 0.0, perror = 0.0;
    double error = 0.0;
    
    
    double x[MAX_NUM_POINTS];
    double y[MAX_NUM_POINTS];
    double corr[MAX_NUM_POINTS];
    
    char datafile[LINE_SIZE];
    sprintf(datafile, "data/response/%s.txt", scan);
    FILE *data = fopen(datafile, "r");
    char correctionfile[LINE_SIZE];
    sprintf(correctionfile, "correction/%s.txt", corr_table);
    FILE *correction = fopen(correctionfile, "r");
    if (!data) {
        printf("Invalid data file.\n");
        return;
    }
    
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    
    if (correction) {
        for (int h = 0; h < corrheader; h++) 
            if (!fgets(line, sizeof(line), correction)) printf("Error");
        while(fgets(line, sizeof(line), correction)) {
            strtok(line, "\n");
            if (index * inc >= num_points) break;
            double xval = strtod(line, &p);
            double yval = atof(p);
            corr[(int) xval - offset] = yval;
            index++;
        }
        for (int i = 0; i < num_points - 1; i += inc) {
            for (int j = i + 1; j < i + inc; j++) {
                double m = (corr[i + inc] - corr[i]) / inc;
                corr[j] = m * (j - i) + corr[i];
            }
        }
    } else if (correctionfile) {
        printf("Correction lookup table not found.\n" );
    }
    
    index = 0;
    for (int h = 0; h < dataheader; h++) 
        if (!fgets(line, sizeof(line), data)) printf("Error");
    while(fgets(line, sizeof(line), data)) {
        strtok(line, "\n");
        if (index >= num_points) break;
        double xval = strtod(line, &p);
        double yval = atof(p);
        x[index] = xval;
        y[index] = yval;
        if ((index < NAP/2) || (index >= num_points - NAP/2)) {
            baseline += yval;
        }
        index++;
    }
    
    fclose(data);
    if (correction) fclose(correction);
    baseline = baseline / (double) NAP;
    for (int j = 0; j < num_points; j++) {
        y[j] -= baseline;
        if (correction) y[j] *= corr[j];
        if ((j + offset > lbmin) && (j + offset < lbmax)) {
            L_b += y[j];
            lberror += pow(corr[j], 2) * (y[j] + baseline + baseline) 
                       + pow(y[j], 2) * pow(.04148 * corr[j], 2);   
        } else if ((j + offset >= pbmin) && (j + offset <= pbmax)) {
            P_b += y[j];
            perror += pow(corr[j], 2) * (y[j] + baseline + baseline) + pow(y[j], 2) * 
                      pow(.04148 * corr[j], 2); 
        }
    }
    
    perror = sqrt(perror);
    lberror = sqrt(lberror);
    
    QY = P_b / (L_a - L_b);
    
    printf("P_b %f \ndivided by \nL_a %f \nminus\nL_b %f: \nQY %f\n", P_b, L_a, L_b, QY);
    
    error = P_b * pow(L_a - L_b, -2) + (L_a + L_b) * P_b * P_b * pow(L_a - L_b, -4);
    error = sqrt(error);
    printf("error %f\n", error);
    
    TGraphErrors* g = new TGraphErrors(num_points, x, y);
    g->SetTitle(title);
    if (correction) g->SetLineColor(2);
    else g->SetLineColor(4);
    g->GetXaxis()->SetTitle(xtitle);
    g->GetYaxis()->SetTitle(ytitle);
    g->GetYaxis()->SetTitleOffset(1.5);
    g->GetXaxis()->SetLimits((double) offset, (double) offset + num_points -1);
    g->Draw();

}

