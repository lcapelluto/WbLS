#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*  Read data from FILE, which has HEADER_SIZE header lines
 *  before data text and NUM_POINTS total data points. Make a
 *  TGraph titled TITLE, x-axis labeled XTITLE, y-axis labeled
 *  YTITLE.
 */
void SimpleGraph(char file[] = "",
                 int header_size = 4, int num_points = 350,
                 char title[] = "Graph",
                 char xtitle[] = "wavelength [nm]", 
                 char ytitle[] = "Intensity [photons/sec]") {
                 
    TGraph* plot = new TGraph();
    FILE *data = fopen(file, "r");
    if (!data) {
        printf("Options: \"filename\", header size, num points, \"title\", \"x-axis label\", \"y-axis label\".\n");
        return;
    }
    
    char line[LINE_SIZE];
    char* p;
    int index = 0;
    for (int h = 0; h < header_size; h++) 
        fgets(line, sizeof(line), data);
    while(fgets(line, sizeof(line), data)) {
        strtok(line, "\n");
        if (index > num_points) break;
        double x = strtod(line, &p);
        double y = atof(p);
        plot->SetPoint(index, x, y);
        index++;
    }
    fclose(data);
    
    plot->SetTitle(title);
    plot->GetXaxis()->SetTitle(xtitle);
    plot->GetYaxis()->SetTitle(ytitle);
    plot->GetYaxis()->SetTitleOffset(1.5);
    plot->SetMarkerColor(4);
    plot->SetMarkerStyle(20);
    plot->SetMarkerSize(.7);
    plot->Draw();
}

/* Max line size for data. */
const int LINE_SIZE = 60;

