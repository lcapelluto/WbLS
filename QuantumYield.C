#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Maximum length of data file line. */
const int LINE_SIZE = 100;
/* Number of header lines in data files. */
const int dataheader = 4;
/* Number of header lines in correction files. */
int corrheader = 0;
/* Increment in wavelengths between given corrections. */
int inc = 0;
/* Wavelength that sample data starts at. */
const int offset_data = 300;
/* Wavelength that L_a data starts at. */
int offset_ref = 320;
/* Number of data points. */
int num_points = 351;
/* 'N'umber of 'A'veraging 'P'oints to establish a baseline. */
const int NAP = 20;
//const int NAP = 12;
/* Maximum number of data points in one scan. */
const int MAX_NUM_POINTS = 360;
/* Array holding correction factors for wavelength index + offset_data. */
double corr[MAX_NUM_POINTS];
/* Type of correction. 0 = no correction, 1 = default, 2 = quanta. */
int is_corr = 0;
/* Important shared values in the quantum yield calculations and their errors.
 *
 * L: Counts of photons in excitation peak of a) solvent without sample
 * b) sample directly in beam c) sample indirectly excited by the beam. */
double L_a = 0.0;
double L_a_err = 0.0;
double L_b = 0.0;
double L_b_err = 0.0;
double L_c = 0.0;
double L_c_err = 0.0;
/* P: Counts of photons in photoluminscence peak with sample b) directly in
 * beam c) not directly in beam. */
double P_b = 0.0;
double P_b_err = 0.0;
double P_c = 0.0;
double P_c_err = 0.0;

// make l range from excitation -8 to excitation +4.
// can also specify file names this way
//solute, solvent, excitation wavelength!
// Automate min and max??

/* Print the quantum yield (QY) values obtained by using the 2MM measurement
 * method and by using the 3MM method. 2MM will be obtained for each of 
 * IN_BEAM_DATA and OUT_DATA that are provided, 3MM will be obtained only if
 * both are provided. LMIN, LMAX, PMIN, and PMAX set bounds on peak integration.
 * Emission correction for the sensitivity of the detector is specified as
 * EMISSION_CORRECTION, and can be either "default" or "quanta" or any prefix
 * of them.
 */
void QuantumYield(char* solvent_data_file = "",
                     char* in_beam_data = "",
                     int lmin = 0, int lmax = 0,
                     int pmin = 0, int pmax = 0,
                     char* out_beam_data = "",
                     char* emission_correction = "d") {
    double QY;
    double error;
    
    FILE* solvent = open_data(solvent_data_file);
    FILE* in_data = open_data(in_beam_data);
    FILE* out_data = open_data(out_beam_data);
    FILE* correction = open_correction(emission_correction);

    if (!solvent || (!in_data && !out_data)) {
        printf("Invalid data file.\n");
        printf("Arguments: reference data file, in-beam data file, lmin, lmax, pmin, pmax, indirect data file, correction name.\n");
        close_files(solvent, in_data, out_data, correction);
        return;
    }

    if (is_corr) {
        read_correction(correction);
    }

    get_L_a(solvent, lmin, lmax);
    
    if (in_data) {
        QY = TwoMMQY(in_data, lmin, lmax, pmin, pmax);
        error = get_error();
        printf("\n\nL_a %f\nL_b %f\nP_b %f\n2MMQY %f +/- %f\n", L_a, L_b, P_b,
                                                              QY, error);
        if (out_data) {
            QY = ThreeMMQY(out_data, lmin, lmax, pmin, pmax);
            error = get_error();
            printf("L_c %f\nP_c %f\n3MMQY %f +/- %f\n", L_c, P_c, QY, error);
        }
    }
    if (out_data) {
        rewind(out_data);
        L_b = 0;
        P_b = 0;
        QY = TwoMMQY(out_data, lmin, lmax, pmin, pmax);
        error = get_error();
        printf("2MMQY (out data) %f +/- %f\n", QY, error);
    }
    printf("\n");
    
    close_files(solvent, in_data, out_data, correction);
}

/* Determine the value for L_a, which is the number of counts in the
 * excitation peak without the sample in the solvent. Data comes from
 * file DATA, points are summed from MIN to MAX.
 */
void get_L_a(FILE* data, int min, int max) {
    double x[MAX_NUM_POINTS];
    double y[MAX_NUM_POINTS];
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    double baseline = 0.0;
    offset_ref = strip_header(data);
    if (!offset_ref) {
        return;
    }
    while(fgets(line, sizeof(line), data)) {
        strtok(line, "\n");
        if (index >= num_points) break;
        double xval = strtod(line, &p);
        double yval = atof(p);
        x[index] = xval;
        y[index] = yval;
        if ((index < (3 * NAP) / 2)||(index >= num_points - (3 * NAP ) / 2)) {
        //if (index >= num_points - 45) {
            baseline += yval;
        }
        index++;
    }
    baseline = baseline / (3 * NAP);
    //baseline /= 45;
    for (int j = 0; j < num_points; j++) {
        y[j] -= baseline;
        if (is_corr) y[j] *= corr[j + offset_ref - offset_data];
        if ((j + offset_ref) > min && (j + offset_ref) < max) {
            L_a += y[j];
            L_a_err += pow(corr[j], 2) * (y[j] + baseline + baseline) 
                       + pow(y[j], 2) * pow(.04148 * corr[j], 2); 
        }
    }
    L_a_err = sqrt(L_a_err);
    
    TGraph* g = new TGraph(80, x, y);
    g->Draw();
}


/* Calculate L_b and P_b for the 2MM QY measurement method. DATA can be a scan
 * of a directly or indirectly excited sample in the integrating sphere. LMIN,
 * LMAX describe the bounds on the summation of the excitation peak. PMIN and
 * PMAX bound the photoluminescence peak summation. Must be called after get_L_a.
 */
double TwoMMQY(FILE* data, int lmin, int lmax, int pmin, int pmax) {
    double x[MAX_NUM_POINTS];
    double y[MAX_NUM_POINTS];
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    double baseline = 0.0;
    
    if (!strip_header(data)) {
        return 0.0;
    }
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
    baseline = baseline / NAP;
    printf("BASELINE %f\n",baseline);
    for (int j = 0; j < num_points; j++) {
        y[j] -= baseline;
        if (is_corr) y[j] *= corr[j];
        if ((j + offset_data > lmin) && (j + offset_data < lmax)) {
            L_b += y[j];
            L_b_err += pow(corr[j], 2) * (y[j] + baseline + baseline) 
                       + pow(y[j], 2) * pow(.04148 * corr[j], 2);   
        } else if ((j + offset_data >= pmin) && (j + offset_data <= pmax)) {
            P_b += y[j];
            P_b_err += pow(corr[j], 2) * (y[j] + baseline + baseline) 
                       + pow(y[j], 2) * pow(.04148 * corr[j], 2); 
        }
        L_b_err = sqrt(L_b_err);
        P_b_err = sqrt(P_b_err);
    }
    
    return P_b / (L_a - L_b);
}


/* Calculate L_c and P_c for the 3MM QY measurement method. Must be called
 * after TwoMMQY has been called on in data. DATA must specify a scan of
 * an indirectly excited sample. LMIN, LMAX, PMIN and PMAX are the same as
 * in TwoMMQY.
 */
double ThreeMMQY(FILE* data, int lmin, int lmax, int pmin, int pmax) {
    double A = 0.0;
    double y[MAX_NUM_POINTS];
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    double baseline = 0.0;
    
    if (!strip_header(data)) {
        return 0.0;
    }
    while(fgets(line, sizeof(line), data)) {
        strtok(line, "\n");
        if (index >= num_points) break;
        double xval = strtod(line, &p);
        double yval = atof(p);
        y[index] = yval;
        if ((index < NAP/2) || (index >= num_points - NAP/2)) {
            baseline += yval;
        }
        index++;
    }
    
    baseline = baseline / NAP;
    for (int j = 0; j < num_points; j++) {
        y[j] -= baseline;
        if (is_corr) y[j] *= corr[j];
        if ((j + offset_data > lmin) && (j + offset_data < lmax)) {
            L_c += y[j];
        } else if ((j + offset_data >= pmin) && (j + offset_data <= pmax)) {
            P_c += y[j];
        }
    }
    
    A = 1 - (L_b / L_c);
    return (P_b - ((1 - A) * P_c)) / (A * L_a);
}


/********************************************************************/
/********************************************************************/
/********************** Helper functions ****************************/
/********************************************************************/
/********************************************************************/

/* Return QY error. */
double get_error() {
    double error;
    if (is_corr != 1) {
        L_a_err = sqrt(L_a);
        L_b_err = sqrt(L_b);
        L_c_err = sqrt(L_c);
        P_b_err = sqrt(P_b);
        P_c_err = sqrt(P_c);
    }
    error = P_b * pow(L_a - L_b, -2) + (L_a + L_b)
                * P_b * P_b * pow(L_a - L_b, -4);
    error = sqrt(error);
    return error;
}

/* Strip off the header from FILE and change NUM_POINTS. Returns 0 if there
 * is an error, offset otherwise. */
int strip_header(FILE* file) {
    char line[LINE_SIZE];
    char *p = NULL;
    for (int h = 0; h < dataheader; h++) {
        if (!fgets(line, sizeof(line), file)) {
            printf("Error reading data file.\n");
            return 0;
        }
        if (h == 1) {
            num_points =  atoi(line);
        } else if (h == 2) {
            p = line;
            strtok(p, ":-");
            p = strtok(NULL, ":-");
        }
    }
    if (!p) {
        return 0;
    }
    return atoi(p);
}

/* Return the baseline from FILE. Fill up array Y with data from FILE. */
double get_base(FILE* file, int[] y, int n_a_p) {
    int index = 0;
    char line[LINE_SIZE];
    while(fgets(line, sizeof(line), file)) {
        strtok(line, "\n");
        if (index >= num_points) break;
        double xval = strtod(line, &p);
        double yval = atof(p);
        y[index] = yval;
        if ((index < NAP/2) || (index >= num_points - NAP/2)) {
            baseline += yval;
        }
        index++;
    }
    
    return baseline / NAP;
}

/* Open data file, given FILENAME in folder data/response/. */
FILE* open_data(char* filename) {
    char datafile[LINE_SIZE];
    sprintf(datafile, "data/%s.txt", filename);
    return fopen(datafile, "r");
}

/* Open correction file CORRECTION either default or quanta. */
FILE* open_correction(char* correction) {
    if (cmpfirst(correction, "default")) {
        is_corr = 1;
        inc = 1;
        return fopen("correction/alignedDefault.txt", "r");
    } else if (cmpfirst(correction, "quanta")) {
        is_corr = 2;
        corrheader = 6;
        inc = 2;
        return fopen("correction/rejected/emcorr-sphere-quanta.txt", "r");
    } else if (cmpfirst(correction, "cyclohexane")
               || cmpfirst(correction, "cyhx")) {
       is_corr = 3;
       corrheader = 0;
       inc = 10;
       return fopen("correction/alignedCyhx.txt", "r");
    } else if (cmpfirst(correction, "normalizer")) {
       is_corr = 4;
       corrheader = 0;
       inc = 1;
       return fopen("correction/normalizer.txt", "r");
    }
    printf("No emission correction is being applied.\n");
    return NULL;
}

/* Read correction data from CORRECTION into correction array. */
void read_correction(FILE* correction) {
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    for (int h = 0; h < corrheader + (offset_data - 300)/inc; h++) {
        if (!fgets(line, sizeof(line), correction)) {
            printf("Error reading correction file.\n");
            return;
        }
    }
    while(fgets(line, sizeof(line), correction)) {
        strtok(line, "\n");
        if (index * inc >= num_points) break;
        double xval = strtod(line, &p);
        double yval = atof(p);
        corr[(int) xval - offset_data] = yval;
        index++;
    }
    for (int i = 0; i < num_points - 1; i += inc) {
        for (int j = i + 1; j < i + inc; j++) {
            double m = (corr[i + inc] - corr[i]) / inc;
            corr[j] = m * (j - i) + corr[i];
        }
    }
}

/* Close 1 - 4 data files F0, F1, F2, F3. */
void close_files(FILE* f0, FILE* f1 = NULL, FILE* f2 = NULL, FILE* f3 = NULL) {
    if (f0) fclose(f0);
    if (f1) fclose(f1);
    if (f2) fclose(f2);
    if (f3) fclose(f3);
}

/* Return 1 if the leading characters in S0 and S1 match until one string
 * ends, 0 otherwise. */
int cmpfirst(char* s0, char* s1) {
    if (!(*s0) || !(*s1)) {
        return 1;
    } else if (*s0 == *s1) {
        return cmpfirst(s0 + 1, s1 + 1);
    } else {
        return 0;
    }
}

