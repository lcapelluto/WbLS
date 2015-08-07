const int MAX_NUM_POINTS = 651;
const int LINE_SIZE = 60;

void Reabsorption() {
    int num_points = 347;
    int NAP = 20;
    double x[MAX_NUM_POINTS];
    double y[MAX_NUM_POINTS];
    int xcorr[MAX_NUM_POINTS];
    double corr[MAX_NUM_POINTS];
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    double xval, yval;
    
    int a = 70;
    int b = 100;
    int c = 135;
    
    FILE *data = fopen("data/reabsorption/anthracene_etoh_30.2_box.txt", "r");
    FILE *correction = fopen("correction/emcorri.txt", "r");
    if (!data) {
        printf("Invalid data file.\n");
        return;
    }

    if (correction) {
        for (int h = 0; h < 31; h++) 
            if (!fgets(line, sizeof(line), correction)) printf("Error");
        while(fgets(line, sizeof(line), correction)) {
            strtok(line, "\n");
            if (index > (350/2)) break;
            xval = strtod(line, &p);
            yval = atof(p);
            index++;
            xcorr[index] = xval;
            corr[(int) xval] = yval;
        }
        for (int i = 1; i < (352/2); i++) {
            double m = (corr[xcorr[i+1]] - corr[xcorr[i]]) /  (xcorr[i+1] - xcorr[i]);
            for (int j = xcorr[i] + 1; j < xcorr[i + 1]; j++) {
                corr[j] = m * (j - xcorr[i]) + corr[xcorr[i]];
            }
        }
        
    }
     
    index = 0;
    double baseline = 0.0;
    for (int h = 0; h < 4; h++) 
        if (!fgets(line, sizeof(line), data)) printf("Error");
    while(fgets(line, sizeof(line), data)) {
        strtok(line, "\n");
        if (index >= num_points) break;
        xval = strtod(line, &p);
        yval = atof(p);
        x[index] = xval;
        y[index] = yval;
        if ((index < NAP/2) || (index >= num_points - NAP/2)) {
            baseline += yval;
        }
        index++;
    }
    
    fclose(data);
    if (correction) fclose(correction);
    
    
    baseline = baseline / NAP;
    double peak1 = 0.0, norm1 = 0.0;
    printf("baseline %f\n", baseline);
    for (int j = 0; j < num_points; j++) {
        y[j] -= baseline;
        if (correction) y[j] *= corr[(int) x[j]];
        if (j > a && j < c) {
            peak1 += y[j];
        }
        if (j > b && j < c) {
            norm1 += y[j];
        }
    }
    
    TGraph *FS = new TGraph(300, x, y);
    
    ///******************************************************************//
    ///******************************************************************//
    ///******************************************************************//
    ///******************************************************************//
    
    FILE *data2 = fopen("data/response/anthracene_ethanol/AE_340_IN_default.txt", "r");
    FILE *correction2 = fopen("correction/alignedCyhx.txt", "r");
    if (!data2) {
        printf("Invalid data file2.\n");
        return;
    }
    index = 0, baseline =0.0;
    double bkgrd[MAX_NUM_POINTS];
    FILE *background = fopen("data/WbLS/DayaBay_250.txt","r");
    if (background) {
        while(fgets(line, sizeof(line), background)) {
            strtok(line, "\n");
            xval = strtod(line, &p);
            yval = atof(p);
            bkgrd[index] = yval;
            index++;
            if ((index < NAP/2) || (index >= num_points - NAP/2)) {
                baseline += yval;
            }
        }
    }
    baseline = baseline / NAP;
    printf("baseline %f\n", baseline);
    index = 0;

    if (correction2) {
        while(fgets(line, sizeof(line), correction2)) {
            strtok(line, "\n");
            xval = strtod(line, &p);
            yval = atof(p);
            index++;
            corr[(int) xval] = yval;
        }        
    }    
    
    double x2[MAX_NUM_POINTS];
    double y2[MAX_NUM_POINTS];
    
    index = 0;
    double baseline = 0.0;
    for (int h = 0; h < 4; h++) 
        if (!fgets(line, sizeof(line), data2)) printf("Error");
    while(fgets(line, sizeof(line), data2)) {
        strtok(line, "\n");
        if (index >= num_points) break;
        xval = strtod(line, &p);
        yval = atof(p);
        x2[index] = xval;
        y2[index] = yval;
        if ((index < NAP/2) || (index >= num_points - NAP/2)) {
            baseline += yval;
        }
        index++;
    }
    
    fclose(data2);
    if (correction2) fclose(correction2);
    baseline = baseline / NAP;
    double peak2 = 0.0, norm2 = 0.0;
    printf("baseline %f\n", baseline);
    for (int j = 0; j < num_points; j++) {
        //y2[j] -= (bkgrd[j]-2900);
        y2[j] -= baseline;
        if (correction2) y2[j] *= corr[(int) x2[j]];
        if (j > (a-3) && j < (c-3)) {
            peak2 += y2[j];
        }
        if (j > (b-3) && j < (c-3)) {
            norm2 += y2[j];
        }
        y2[j] *= 7500;
    }

    TGraph *IS = new TGraph(300, x2, y2);
    IS->Draw();
    FS->Draw("SAME");
    
    double factor = norm1/norm2;
    printf("norm1 %f, norm2 %f, factor %f \n", norm1, norm2, factor);
    printf("peak2 %f / peak1 %f = %f\n", peak2 * factor, peak1, (peak2*factor)/ peak1);
    
}

