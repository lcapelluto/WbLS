#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Array size for x and y points. */
const int MAX_NUM_POINTS = 500;
/* Estimated maximum length of data file line. */
const int LINE_SIZE = 60;

void Fill() {
    int xcorr[400];
    double corr[650];
    
    FILE* r = fopen("CyhxCorrectionDataPeak.txt", "r");
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    if (r) {
        while(fgets(line, sizeof(line), r)) {
            strtok(line, "\n");
            double xval = strtod(line, &p);
            double yval = atof(p);
            xcorr[index] = (int) xval;
            corr[(int) xval] = yval;
            index++;
        }
        for (int i = 0; i < 37; i++) {
            double m = (corr[xcorr[i+1]] - corr[xcorr[i]]) /  (xcorr[i+1] - xcorr[i]);
            for (int j = xcorr[i] + 1; j < xcorr[i + 1]; j++) {
                printf("i %d j %d y %f\n", i, j, corr[xcorr[i]]);
                corr[j] = m * (j - xcorr[i]) + corr[xcorr[i]];
            }
        }
        
    }
    
    FILE* w = fopen("alignedCyhx.txt","w");
    for (int i = 300; i < 648; i++) {
        printf("%dnm - %f\n", i, corr[i]);
        sprintf(line, "%d %f\n", i, corr[i]);
        fputs(line, w);
    } 

    fclose(r);
    fclose(w);
}

