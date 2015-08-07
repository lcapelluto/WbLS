#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Array size. */
const int MAX_NUM_POINTS = 650;
/* Estimated maximum length of data file line. */
const int LINE_SIZE = 60;
/* Number of Averaging Points to establish a baseline. */
const int NAP = 12;

/* Useful for corrections with inconsistent intervals in wavelength.
 * Apply default emission correction, subtract baseline, and make 
 * a TGraph with NUM_POINTS points of the resulting response from DATAFILE,
 * titled TITLE with x and y axes labeled XTITLE and YTITLE, respectively. Data
 * must have DATAHEADER number of header lines. If LOW? is set to 1, it will
 * assume excitation peak is too close to 300nm to use bottom values for
 * baseline.
 */
void GraphResponse(char* datafile = "", char* title = "Graph",
              int num_points = 31, int dataheader = 4,
              int is_low = 0,
              char* xtitle = "wavelength [nm]",
              char* ytitle = "Intensity [photons/sec]") {


    double x[MAX_NUM_POINTS];
    double y[MAX_NUM_POINTS];
    int xcorr[MAX_NUM_POINTS];
    double corr[MAX_NUM_POINTS];
    
    FILE *data = fopen(datafile, "r");
    FILE *correction = fopen("correction/water.txt", "r");
    if (!data) {
        printf("Invalid data file.\n");
        return;
    }
    
    char line[LINE_SIZE];
    char* p;
    int index = 0;

    if (correction) {
        while(fgets(line, sizeof(line), correction)) {
            strtok(line, "\n");
            double xval = strtod(line, &p);
            double yval = atof(p);
            index++;
            xcorr[index] = xval;
            corr[(int) xval] = yval;
        }
        for (int i = 1; i < 38; i++) {
            double m = (corr[xcorr[i+1]] - corr[xcorr[i]]) /  (xcorr[i+1] - xcorr[i]);
            for (int j = xcorr[i] + 1; j < xcorr[i + 1]; j++) {
                corr[j] = m * (j - xcorr[i]) + corr[xcorr[i]];
            }
        }
        
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
        printf("x %f, y %f\n",xval,yval);
        if (((index < NAP/2) || (index >= num_points - NAP/2)) && !is_low) {
            baseline += yval;
        } else if ((index < NAP) && (is_low)) {
            baseline += yval;
        }
        index++;
    }
    
    fclose(data);
    if (correction) fclose(correction);
    baseline = baseline / NAP;
    printf("baseline %f\n", baseline);
    for (int j = 0; j < num_points; j++) {
        y[j] -= baseline;
        if (correction) y[j] *= corr[(int) x[j] + 2];
    }
    
    
    
    TGraphErrors* g = new TGraphErrors(num_points, x, y);
    g->SetTitle(title);
    g->GetXaxis()->SetTitle(xtitle);
    g->GetYaxis()->SetTitle(ytitle);
    g->GetYaxis()->SetTitleOffset(1.5);
    g->Draw();

}

