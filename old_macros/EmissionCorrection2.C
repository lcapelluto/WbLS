#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Assign the correction factor to the proper wavelength, accounting for
 * instrumentation offsets.
 *
 * Determine the correction to be applied to emission spectra from the
 * integrating sphere, by printing 3 arrays holding the inverse of the sum
 * (with an averaged baseline subtracted from each point) of events from each
 * scan.
 *
 * Data was collected using a high pressure xenon lamp with two monochromators 
 * to select one exitation wavelength, which traveled through a light tube
 * to an integrating sphere, empty except for the baffle and sample holder.
 * Reflected light went through another light tube and was recorded by a 
 * third monochromator with FeliXGX data acquisition software.
 * 
 * 73 emission scans were taken, with a +-15nm range from the excitation 
 * wavelength. Gain was set to 6.8V, slit widths were all 1mm, integration
 * time was 10 seconds. Each run has uncorrected raw data and an excitation
 * correction given by FeliX, either the Default Correction or excorr.
 */
void EmissionCorrection2() {
    double default_corr[NUM_CORR];
    double dc_x[NUM_CORR];

    
    char filename[LINE_SIZE];    
    char id[] = "02";
    char wavelength[] = "300";
    char correction[] = "default";    
    for (int i = 0; i < NUM_FILES; i++) {
        double offset = 0.0;
        double arr_y[NUM_POINTS + 1];
        sprintf(filename, "data/empty-sphere/%s_%snm_%s.txt", id, wavelength, 
                                                         correction);
        FILE *data = fopen(filename, "r");
        
        char line[LINE_SIZE];
        char* p;
        int index = 0;
        double sum = 0.0;
        
        double max = 0;
        double peak = 0;
        for (int h = 0; h < HEADER_SIZE; h++) 
            fgets(line, sizeof(line), data);
            
        while(fgets(line, sizeof(line), data)) {
            strtok(line, "\n");
            if (index >= NUM_POINTS) break;
            double x = strtod(line, &p);
            double y = atof(p);
            if ((index < NUM_AVG) || (index >= NUM_POINTS - NUM_AVG)) {
                offset += y;
            }
            if (y > max) {
                peak = x;
                max = y;
            }
            arr_y[index] = y;
            index++;

        }
        fclose(data);
        
        offset /= (NUM_AVG * 2);
        for (int j = 0; j < NUM_POINTS; j++) {
            arr_y[j] -= offset;
            sum += arr_y[j];
        }
        /* Tried to normalize around 430nm. */
        if (!strcmp(correction, "default")) {
            default_corr[i/4] = 1/(sum * 0.0000007);
            dc_x[i/4] = peak;
            printf("peak %d\ncorr %f\n",peak,default_corr[i/4]);
        }
        
        updateFilename(i, id, wavelength, correction); 
        
    }
    
    FILE* default_result = fopen("correction/defaultCorrectionData2.txt", "w");

    printf("\n__________________\nDEFAULT 1/sum \n");
    for (int k = 0; k < NUM_CORR - 1; k++) {
        printf("%dnm: %f\n", dc_x[k], default_corr[k]);
        sprintf(line, "%d %f\n", dc_x[k], default_corr[k]);
        fputs(line, default_result);
    }

    fclose(default_result);
    
    /* VARIATION.C was used to calculate the average fraction of the
     * standard error over the mean photon count for each of the three
     * excitation corrections. This is the factor used for error bars.
     */

    double def_err[NUM_CORR];

    for (int n = 0; n < NUM_CORR - 1; n++) {
        def_err[n] = .01371844 * default_corr[n];
    }

    TGraphErrors* defaultG = new TGraphErrors(36, dc_x, default_corr, NULL, def_err);
    defaultG->SetLineColor(7); //cyan
    defaultG->Draw("SAME");

}

/* Modifies ID, WAVELENGTH and CORRECTION strings to open the next data file
 * using index parameter I. */
void updateFilename(int i, char* id, char* wavelength, char* correction) {
    if ((i + 1) % 4 == 0) sprintf(wavelength, "%d", atoi(wavelength) + 10);
    if (i % 2 != 0) sprintf(id, "%02d", atoi(id) + 1);
    if (i % 2 == 0) {
        strcpy(correction, "uncorr");
    } else if ((((int) atoi(id) + 1) % 2) != 0) {
        strcpy(correction,"default");
    } else {
        strcpy(correction, "excorr");
    }
    
}

/* Constants regarding data acquisition. */
const int NUM_POINTS = 31;
const int NUM_FILES = 144;
/* Number of header lines in FeliXGX generated data text files. */
const int HEADER_SIZE = 4;
/* Constants regarding array sizes. */
const int LINE_SIZE = 60;
const int NUM_CORR = 37;
const int NUM_UNCORR = 73;
/* Number of points used on each side of curve to determine offset. */
const int NUM_AVG = 3;

