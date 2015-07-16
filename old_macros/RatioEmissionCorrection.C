#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Constants regarding data acquisition. */
const int NUM_POINTS = 151;
/* Constants regarding array sizes. */
const int LINE_SIZE = 60;

/* Determine the correction to be applied to emission spectra by taking a ratio
 * of an accepted reference to data that was measured in the integrating sphere.
 */
void RatioEmissionCorrection() {
    double y[NUM_POINTS + 1];
    double corr[NUM_POINTS];

    FILE *data = fopen("data/ratio emcorr/POPOP_IS.txt", "r");
    FILE *ref = fopen("data/ratio emcorr/POPOPref.txt", "r");
    
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    while(fgets(line, sizeof(line), data)) {
        strtok(line, "\n");
        double xval = strtod(line, &p);
        double yval = atof(p);
        y[index] = yval - 3424.978;
        index++;
    }
    
    index = 0;
    while(fgets(line, sizeof(line), ref)) {
        strtok(line, "\n");
        double xval = strtod(line, &p);
        double yval = atof(p);
        corr[index] = (yval - 3500) / y[index];
        index++;
        fgets(line, sizeof(line), ref);
        strtok(line, "\n");
    }
    
    FILE* result = fopen("correction/ratioCorrectionData.txt", "w");
    for (int k = 0; k < NUM_POINTS; k++) {
        printf("%dnm: %f\n", k + 350, corr[k]);
        sprintf(line, "%d %f\n", k + 350, corr[k]);
        fputs(line, result);
    }
    
    fclose(data);
    fclose(ref);
    fclose(result);
    
    /* VARIATION.C was used to calculate the average fraction of the
     * standard error over the mean photon count for each of the three
     * excitation corrections. This is the factor used for error bars.
     */
    double err[NUM_POINTS + 1];
    double w[NUM_POINTS + 1];
    double y_corr[NUM_POINTS + 1];
    for (int n = 0; n < NUM_POINTS; n++) {
        err[n] = 0.04148067 * corr[n];
        w[n] = n + 350;
        y_corr[n] = y[n] * corr[n];
    }
    
    
    TGraphErrors* g = new TGraphErrors(151, w, corr, NULL, err);
    g->SetTitle("Emission Corrections");
    g->SetLineColor(4);
    g->Draw();
} 

