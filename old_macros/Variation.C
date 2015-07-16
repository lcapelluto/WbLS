#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "TCanvas.h"

/* Number of data points. */
const int NUM_POINTS = 31;
/* Range of emission scan. */
const double RANGE[] = {415.0, 445.0};
/* Number of header lines in data text files. */
const int HEADER_SIZE = 4;
/* Constants regarding array sizes. */
const int LINE_SIZE = 60;

/* Create a TGraphErrors in analysis of the deviation between scans,
 * using standard error. Print std error, std deviation for each point.
 * Excitation CORRECTION can be "default", "excorr", or "uncorr", the
 * first two being lookup tables supplied by FeliXGX software, the last
 * being uncorrected, raw data.
 *
 * All scans took place one after the other, with excitation of 430nm,
 * 1mm slit width, emission range of 415-445 nm, 6.8V gain, 10 second
 * integration time.
 */
void Variation(char* correction = "default") {
    TGraphErrors* g = new TGraphErrors(NUM_POINTS);
    double average = 0.0;
    for (int w = RANGE[0]; w <= RANGE[1]; w++) {
        average += compute_percent(w, g, correction);
    }
    g->SetTitle("Standard deviation of the excitation peak");
    g->GetXaxis()->SetTitle("wavelength [nm]");
    g->GetYaxis()->SetTitle("Intensity [photons/sec]");
    g->GetXaxis()->SetLimits(RANGE[0], RANGE[1]);
    g->Draw();
    
    printf("Average stdev/mean * 100 = %f\n", average/NUM_POINTS);
}

/* Return the standard deviation of data taken at WAVELENGTH as a
 * percent of total average count at that wavelength. The average
 * count at that wavelength is plotted on PLOT with error given by
 * the standard error.
 */
double compute_percent(int wavelength, TGraphErrors* plot, char correction[]) {
    char filename[LINE_SIZE];
    char id[3];
    int num_files;
    if (!strcmp(correction, "default")) {
        strcpy(id, "01");
        num_files = 8;
    } else if (!strcmp(correction, "uncorr")) {
        strcpy(id, "01");
        num_files = 14;
    } else {
        strcpy(id, "07");
        num_files = 6;
    }
    double sum = 0.0;
    double vals[15];
    
    for (int i = 0; i < num_files; i++) {
        sprintf(filename, "data/variation/%s_%s.txt", id, correction);
        FILE *data = fopen(filename, "r");
        
        char line[LINE_SIZE];
        char* p;
        int index = 0;
        for (int h = 0; h < HEADER_SIZE; h++) 
            fgets(line, sizeof(line), data);
        while(fgets(line, sizeof(line), data)) {
            strtok(line, "\n");
            if (index >= NUM_POINTS) break;
            long int x = strtol(line, &p, 10);
            double y = atof(p);
            if (x == wavelength) {
                sum += y;
                vals[i] = y;
            }
            index++;
        }
        fclose(data);
        
        /* Update file name. */
        sprintf(id, "%02d", atoi(id) + 1);
    }
    
    double mean = sum / num_files;
    double std_error = 0.0;
    for (int i = 0; i < num_files; i++) {
        std_error += pow(vals[i] - mean, 2);
    }
    std_error /= (num_files - 1);
    double stdev = sqrt(std_error);
    std_error = stdev / sqrt(num_files);
    
    printf("%d nm : sd/mean %f se/mean %f\n", wavelength, stdev/mean *100,
                                              std_error/mean * 100);
    plot->SetPoint(wavelength - RANGE[0], wavelength, mean);
    plot->SetPointError(wavelength - RANGE[0], 0.0, std_error);
    
    return stdev/mean *100;
}

