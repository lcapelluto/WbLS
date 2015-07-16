#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Constants regarding data acquisition. */
const int NUM_POINTS = 36;
/* Constants regarding array sizes. */
const int LINE_SIZE = 60;

/* Compare cyhx emcorr with default emcorr.
 */
void ratio() {
    double x[NUM_POINTS + 1];
    double y[NUM_POINTS + 1];
    double corr[NUM_POINTS];

    FILE *data = fopen("correction/CyhxCorrectionData2.txt", "r");
    FILE *ref = fopen("correction/defaultCorrectionData.txt", "r");
    
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    while(fgets(line, sizeof(line), data)) {
        strtok(line, "\n");
        double xval = strtod(line, &p);
        double yval = atof(p);
        x[index] = xval;
        y[index] = yval;
        index++;
    }
    
    index = 0;
    while(fgets(line, sizeof(line), ref)) {
        strtok(line, "\n");
        double xval = strtod(line, &p);
        double yval = atof(p);
        corr[index] = y[index]/yval;
        index++;
    }

    /* VARIATION.C was used to calculate the average fraction of the
     * standard error over the mean photon count for each of the three
     * excitation corrections. This is the factor used for error bars.
     */
    double err[NUM_POINTS + 1];
    double y_corr[NUM_POINTS + 1];
    for (int n = 0; n < NUM_POINTS; n++) {
        err[n] = 0.04148067 * corr[n];
        y_corr[n] = y[n] * corr[n];
    }
    
    
    TGraphErrors* g = new TGraphErrors(36, x, corr, NULL, err);
    g->SetTitle("Ratio of Cyhx Emcorr to Default Emcorr");
    g->SetLineColor(4);
    g->Draw();
} 

