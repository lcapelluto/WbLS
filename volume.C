#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Creates a plot of volume against integrated intensity.
 * 
 * The quartz vial was filled up with DI water then scanned
 * twice with a 30 nm range, 10 sec integration, and excitation
 * 350 then 320. Water was repeatedly removed from the vial and
 * then the scans were performed again.
 *
 * The purpose is to understand the importance of filling the
 * cuvette with the same volume of solvent only and solvent with
 * the fluorophore. I would have expected that removing water would
 * increase the integrated peak intensity, but it seems that there
 * is no such trend.
 */
void volume() {    
    char filename[LINE_SIZE];
    char line[LINE_SIZE];
    double volumes[] = {7.64740, 8.41104, 8.86870, 9.11126, 9.26344, 9.36596};  
    TGraph *g = new TGraph(NUM_FILES);
   
    for (int i = 0; i < NUM_FILES; i++) {
        double offset = 0.0;
        double y[NUM_POINTS + 1];
        int index = 0;
        char* p;
        double sum = 0.0;

        sprintf(filename, "data/vol/%.5f_320.txt", volumes[i]);
        printf("%s\n\n",filename);
        FILE *data = fopen(filename, "r");
        for (int h = 0; h < HEADER_SIZE; h++) 
            fgets(line, sizeof(line), data);      
        while(fgets(line, sizeof(line), data)) {
            strtok(line, "\n");
            if (index >= NUM_POINTS) break;
            double xval = strtod(line, &p);
            double yval = atof(p);
            if ((index < NUM_AVG) || (index >= NUM_POINTS - NUM_AVG)) {
                offset += yval;
            }
            y[index] = yval;
            index++;
        }
        fclose(data);
        
        offset /= (NUM_AVG * 2);
        for (int j = 0; j < NUM_POINTS; j++) {
            y[j] -= offset;
            if (((j+305) >= 313) && ((j+305) <=325)) {
                sum += y[j];
            }
        }
        g->SetPoint(i, volumes[i] - 6.28759, sum);
    }
    g->Draw();

}

/* Constants regarding data acquisition. */
const int NUM_POINTS = 31;
const int NUM_FILES = 6;
/* Number of header lines in FeliXGX generated data text files. */
const int HEADER_SIZE = 4;
/* Constants regarding array sizes. */
const int LINE_SIZE = 60;
/* Number of points used on each side of curve to determine offset. */
const int NUM_AVG = 6;

