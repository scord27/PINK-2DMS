# Brief user guide
Needed packages can be found in requirements.txt and installed using the install.bat file. GUIs for each section can be accessed either through python scripts or through the batch files.
## P.I.N.K Pick
1. Import data "Import surface plot", can change between .metal and .msh5 through dropdown menu. 2DMS data processed through SPIKE requires the method file to also be selected, and if recalibrated ML1-3 values are known for SPIKE processed files these should be type in prior to importing the file.
2. Set mass limits and low threshold for peak picking.
3. Find local maxima across the surface plot between the mass limits, adjust mass limits and low threshold if needed and refind local maxima.
4. Set peak picking grid size and perform peak picking using a 2D Gaussian model (slices each 3D peak through the local maxima and models with a Gaussian) or 3D Gaussian model (model full 3D peak using a 2D Gaussian)
5. Save 3D peak list, change the result to save by dropdown menu (either local maxima, 2D Gaussian, or 3D Gaussian)
## P.I.N.K Sort
1. Import 3D peak list as .csv file, "Pick your .csv".
2. Set sorting parameters at a reference *m/z* value, this will be adjusted across the Y-axis mass range for dynamic binning. One standard deviation is typically a suitable jitter value.
3. Sort the 3D peaklist into MS<sup>2</sup> peak lists, "Sort". Identified precursors will be displayed on initial screen.
4. MS<sup>2</sup> are viewed on the MSMS spectrum tab by double clicking on an entry from the list of precursors.
5. NL and CFLs can be displayed in the "Neutral loss and common fragment spectrum" tab by entering the value, selecting the type by dropdown menu and pressing "PLOT".
6. The full list of MS<sup>2</sup> spectra can be exported as a .mgf or .mzml file on the "exporting as..." tab, chnage between formats using the dropdown menu and press "export as..."
## P.I.N.K Match
1. Import mzML file of MS<sup>2</sup> spectra, "Import mzML files". The MS2 will be displayed (double click) and can be removed (double right click) by selecting from the precursor *m/z* list
2. Import database files, .msp files are compatible.
3. Set ppm error for matching peaks for spectral matching and search databases.
4. Potential assignments can be viewed the "View assignments" tab by double clicking the precursor *m/z* list and double clicking potential matches (which are scored between -1 and 1).
5. Confirm assignments by double right clicking which will be stored in the "View results" tab and saved as a .csv.
## P.I.N.K PCA
1. Import multiple .csv files of 2DMS 3D peak lists (minimum of 3), "Pick your peaklist", ensuring each has a unique name. If peaks need to be filtered out (e.g. a contaminant) this an be done by adding the *m/z* value to the filterpeaks.csv file in resources. Double clicking on a dataset displays the scatter plot double right clicking removes the dataset.
2. Set jitter in the X- and Y-axis at a reference *m/z* value, these are used to match peaks across datasets and are dynamic with *m/z*.
3. Perform PCA analysis, results are shown on the second tab as a 2D score plot, 3D score plot, scree plot and loadings plot (for the 2D score plot).
## Parameter estimator
1. Enter key parameters for your instrument, magnetic field strength in T and low excitation frequency (can be found in method files under EXC_Freq_Low).
2. Enter lowest precursor *m/z* of interest and press "Go". The highest time delay increment you can use for 2DMS will be returned.
3. Enter selected time delay increment and magnetic field strength in T on second tab.
4. Enter a pair of precursor ions *m/z* values and press "Go". A predicted required number of scans to separate the pair by at least one standard deviation will be returned. This can be repeated for multiple pairs if needed.
