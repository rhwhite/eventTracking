# eventTracking
Scripts to process and analyse event-based precipitation data created by the Forward-in-time (FiT) tracking algorithm (Skok et al. 2009). These scripts were used to produce the results published in White et al. 2017.

These scripts are work in progress. To my knowledge they should all work given the right input files, and output directories. The input files can be found on my webpage at http://atmos.washington.edu/~rachel/precipevents.html . You will need to update directories based on your file system. Feel free to copy parts or all of the code, but you should implement your own checks to ensure the code works as you expect. Please let me know of any bugs you find.

If you wish to run the FiT algorithm yourself, the code can be obtained by contacting Gregor Skok at the University of Ljubljana. Depending on your data, you may need to preprocess it. You will then need to use the thresholding script in the preprocess folder to create input files for the algorithm. You will need to create a valid namelist, with filename equal to that expected in run_all.sh. run_all.sh should run the algorithm and perform post-processing, before producing a couple of basic analysis plots. 

Skok et al. 2009. Object-Based Analysis of Satellite-Derived Precipitation Systems over the Low- and Midlatitude Pacific Ocean, Monthly Weather Review, 137(10), 3196–3218, doi:10.1175/2009MWR2900.1
White et al. 2017. Tracking precipitation events in time and space in gridded observational data, in press, GRL
