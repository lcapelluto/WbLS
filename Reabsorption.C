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
    
    FILE *data = fopen("data/reabsorption/PPO_etoh_3.15_box.txt", "r");
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
        if (j > 43 && j < 142) {
            peak1 += y[j];
        }
        if (j > 77 && j < 142) {
            norm1 += y[j];
        }
    }
    
    
    ///******************************************************************//
    ///******************************************************************//
    ///******************************************************************//
    ///******************************************************************//
    
    FILE *data2 = fopen("data/reabsorption/PPO_etoh_301_IS.txt", "r");
    FILE *correction2 = fopen("correction/alignedDefault.txt", "r");
    if (!data2) {
        printf("Invalid data file2.\n");
        return;
    }

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
    
    index = 0;
    double baseline = 0.0;
    for (int h = 0; h < 4; h++) 
        if (!fgets(line, sizeof(line), data2)) printf("Error");
    while(fgets(line, sizeof(line), data2)) {
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
    
    fclose(data2);
    if (correction2) fclose(correction2);
    baseline = baseline / NAP;
    double peak2 = 0.0, norm2 = 0.0;
    printf("baseline %f\n", baseline);
    for (int j = 0; j < num_points; j++) {
        y[j] -= baseline;
        if (correction2) y[j] *= corr[(int) x[j]];
        if (j > 41 && j < 140) {
            peak2 += y[j];
        }
        if (j > 75 && j < 140) {
            norm2 += y[j];
        }
    }
    
    double factor = norm1/norm2;
    printf("factor %f \n", factor);
    printf("peak2 %f / peak1 %f = %f\n", peak2 * factor, peak1, (peak2*factor)/ peak1);
    
}

