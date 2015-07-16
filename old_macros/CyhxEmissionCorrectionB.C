#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* This correction sets the max of each peak to be approximately equal to the
 * others. It is meant as a small experiment in response to an unfortunate
 * pattern noticed in the scans after the correction is applied. 
 *
 * Determine the correction to be applied to emission spectra from the
 * integrating sphere, by printing the inverse of the sum (with an averaged 
 * baseline subtracted from each point) of events from each scan.
 *
 * Data was collected using a high pressure xenon lamp with two monochromators 
 * to select one exitation wavelength, which traveled through a light tube
 * to an integrating sphere which contained a borosilicate vial of cyclohexane.
 * Reflected light went through another light tube and was recorded by a PMT 
 * behind an emission monochromator with FeliXGX data acquisition software.
 * 
 * 39 emission scans were taken, with a +-15nm range from the excitation 
 * wavelength. Gain was set to 6.8V, slit widths were all 1mm, integration
 * time was 10 seconds. Default excitation correction was previously applied
 * to data.
 *
 * An offset between reported wavelengths has been observed: the sample is 
 * excited at the desired excitation wavelength minus ~3 nm, which the emission
 * monochromator reports as -2 nm.
 */
void CyhxEmissionCorrectionB() {
    double corr[NUM_FILES];
    double w[NUM_FILES];
  
    char filename[LINE_SIZE];    
    int wavelength = 300;

    for (int i = 0; i < NUM_FILES; i++) {    
        double baseline = 0.0;
        double x[NUM_POINTS + 1];
        double y[NUM_POINTS + 1];
        sprintf(filename, "data/emcorr w cyhx/%d.txt", wavelength);
        FILE *data = fopen(filename, "r");
        
        char line[LINE_SIZE];
        char* p;
        int index = 0;
        double xmax = 0, ymax = 0;
        for (int h = 0; h < HEADER_SIZE; h++) 
            fgets(line, sizeof(line), data);    
        while(fgets(line, sizeof(line), data)) {
            strtok(line, "\n");
            if (index >= NUM_POINTS) break;
            double xval = strtod(line, &p);
            double yval = atof(p);
            x[index] = xval;
            y[index] = yval;
            if (yval > ymax) {
                xmax = xval;
                ymax = yval;
            }

            if (index < NUM_AVG || index >= NUM_POINTS - NUM_AVG) {
                baseline += yval;
            }
            index++;

        }
        fclose(data);
        baseline /= (NUM_AVG * 2);
        corr[i] = 250000 / (ymax - baseline);
        printf("xmax %d baseline %f\n", xmax, baseline);
        w[i] = xmax;
        
        updateFilename(&wavelength);
    }
    
    FILE* result = fopen("correction/CyhxCorrectionDataB.txt", "w");
    for (int k = 0; k < NUM_FILES; k++) {
        printf("%dnm: %f\n", w[k], corr[k]);
        sprintf(line, "%d %f\n", w[k], corr[k]);
        fputs(line, result);
    }
    fclose(result);
    
    /* VARIATION.C was used to calculate the average fraction of the
     * standard deviation over the mean photon count for each of the three
     * excitation corrections. This is the factor used for error bars.
     */
    double err[NUM_FILES + 1];

    for (int n = 0; n < NUM_FILES; n++) {
        err[n] = 0.04148067 * corr[n];
    }
    
    TGraphErrors* g = new TGraphErrors(39, w, corr, NULL, err);
    g->SetTitle("Emission Corrections");
    g->SetLineColor(4);
    g->Draw();
   
} 

/* Modifies WAVELENGTH to open the next data file. */
void updateFilename(int* wavelength) {
    if (*wavelength < 330) {
        *wavelength += 5;
    } else {
        *wavelength += 10;
    }
}

/* Constants regarding data acquisition. */
const int NUM_POINTS = 31;
const int NUM_FILES = 39;
/* Number of header lines in FeliXGX generated data text files. */
const int HEADER_SIZE = 4;
/* Constants regarding array sizes. */
const int LINE_SIZE = 60;
/* Number of points used on each side of curve to determine baseline. */
const int NUM_AVG = 6;

