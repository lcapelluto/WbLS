#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Array size for x and y points. */
const int MAX_NUM_POINTS = 500;
/* Estimated maximum length of data file line. */
const int LINE_SIZE = 60;
/* Number of Averaging Points to establish a baseline. */
const int NAP = 16;

/* Apply emission correction from CORRECTIONFILE, subtract baseline, and make 
 * a TGraph with NUM_POINTS points of the resulting response from DATAFILE,
 * titled TITLE with x and y axes labeled XTITLE and YTITLE, respectively. Data
 * must have DATAHEADER number of header lines. The correction file has
 * CORRHEADER header lines and values are given in increments of INC.
 */
void Response(char* datafile = "", int offset = 300, char* title = "Graph",
              int num_points = 31, int dataheader = 4,
              char* correctionfile = "correction/alignedCyhx.txt",
              int corrheader = 0, int inc = 1,
              char* xtitle = "wavelength [nm]",
              char* ytitle = "Intensity [photons/sec]") {


    double x[MAX_NUM_POINTS];
    double y[MAX_NUM_POINTS];
    double corr[MAX_NUM_POINTS];
    
    FILE *data = fopen(datafile, "r");
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
            if (index * inc > 350) break;
            double xval = strtod(line, &p);
            double yval = atof(p);
            corr[(int) xval - 300] = yval;
            index++;
        }
        for (int i = 0; i < 350; i += inc) {
            for (int j = i + 1; j < i + inc; j++) {
                double m = (corr[i + inc] - corr[i]) / inc;
                corr[j] = m * (j - i) + corr[i];
                printf("j %d, corr[j] %f\n", j, corr[j]);
            }
        }
    } else if (correctionfile) {
        printf("Correction lookup table not found.\n" );
    }
    
    index = 0;
    double baseline = 0.0;
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
        printf("x %d y %f\n", xval, yval);
    }
    
    fclose(data);
    if (correction) fclose(correction);
    baseline = baseline / NAP;
    printf("baseline %f\n", baseline);
    for (int j = 0; j < num_points; j++) {
        y[j] -= baseline;
        if (correction) y[j] *= corr[(int) x[j] - 300];
        printf("x %d y %f,%d %f\n", x[j] - 300,y[j],j, corr[(int) x[j] - 300]);
    }
    
    
    
    TGraphErrors* g = new TGraphErrors(num_points, x, y);
    g->SetTitle(title);
    if (correction) g->SetLineColor(2);
    else g->SetLineColor(4);
    g->GetXaxis()->SetTitle(xtitle);
    g->GetYaxis()->SetTitle(ytitle);
    g->GetYaxis()->SetTitleOffset(1.5);
    g->Draw();

}

