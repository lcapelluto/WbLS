#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Array size for x and y points. */
const int MAX_NUM_POINTS = 360;
/* Estimated maximum length of data file line. */
const int LINE_SIZE = 60;
/* 'N'umber of 'A'veraging 'P'oints to establish a baseline. */
const int NAP = 60;

/* Graph the emission spectra from data files, apply the selected
 * emission correction, subtract baseline and take a sum, which is L_a.
 * Print L_a. 
 */
void L_a(char* scan = "ethanol/ethanol_340_default",
              int min = 333, int max = 344,
              char* corr_table = "",
              int corrheader = 0, int inc = 10, int offset = 300,
              int num_points = 81, int dataheader = 4) {

    double L_a = 0.0;
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
        printf("Invalid data file.\nUsage: \"datafile\", min, max, \"correctionfile\", corrheader, inc, offset, num points, dataheader.\n");
        return;
    }
    
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    
    if (correction) {
        for (int h = 0; h < corrheader + (offset - 300)/inc; h++) 
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
    } else {
        printf("Correction lookup table not found.\n" );
    }
    
    index = 0;
    double baseline = 0.0;
    for (int h = 0; h < dataheader; h++) 
        if (!fgets(line, sizeof(line), data)) printf("Error\n");
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
        if ((j + offset) > min && (j + offset) < max) {
            L_a += y[j];
            printf("w %d y %f corr %f\n", j + offset, y[j], corr[j]);
            error += pow(corr[j], 2) * (y[j] + baseline + baseline) + pow(y[j], 2) * 
                      pow(.04148 * corr[j], 2); 
        }
    }
    TGraphErrors* g = new TGraphErrors(num_points, x, y);
    g->Draw();
    printf("\nL_a of %s using range %d-%d and %s:\n %f photons.\n\n", datafile, 
                                                min, max, correctionfile, L_a);
    printf("Error: %f\n", sqrt(error));

}

