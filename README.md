# RFCC
Rough-Fuzzy Circular Clustering

### Article: 
P. Maji and S. Mahapatra, "Circular Clustering in Fuzzy Approximation Spaces for Color Normalization of Histological Images," in *IEEE Transactions on Medical Imaging*, pp. 1--11, 2020.
doi: 10.1109/TMI.2019.2956944

URL: https://ieeexplore.ieee.org/document/8918273

The codes are written in C. The details regarding the command line variables are as follows:

f: Name of the file containing the image file names
c: Number of clusters
p: Name of the file containing parameter values
m: fuzzifier value
s: size of the rectangle
a: alpha parameter value
r: Name of the reference/template image

COMMAND: 
   gcc RFCC_vM.c -lm -o rfcc
   
   ./rfcc -f ImageNames.txt -c 3 -p parameters.txt -m 2.0 -s 3 -a 0.1 -r Template.ppm
