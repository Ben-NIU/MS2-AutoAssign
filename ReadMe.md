## MS/MS Spectrum Analysis & Auto-Assignment  
Mass spectrometry takes advantage of MS/MS spectra to provide sequence information of peptide. One has to, however, scrutinize each MS/MS spectrum when it comes to locate the residue that's been covalently labeled. The ***motivation*** of writing this app is to achieve automatic MS/MS spectrum assignments and annotations.  

### **in silico** Fragments  
 * Generate a list of b-/y-ion series of a given peptide 
 * Incorporate residue-level modifications 
 * Have the option to include user-specified neutral loss  
 
 **Shown below is an example of the daughter-ion list.**  
 ![Alt text](pics/silico.png?raw=true "Optional Title")  
 
### Peaks Matching    
 * Read mzXML file ONLY  
 * Display MS/MS spectrum and stats (RT, precursor m/z, charge state) with given scan i.d. 

  **Shown below is an example of the stats box.**  
 ![Alt text](pics/stats.png?raw=true "Optional Title")  
 
 * Annotate spectrum with user-defined parameters  
 
  **Shown below is an example of the MS/MS spectrum.**  
 ![Alt text](pics/spectrum.png?raw=true "Optional Title")  
  
 * Display annotated peaks in Table, which allows selective removal of annotations  
 
   **Shown below is an example of the Table that contains all the annotated peaks.**  
 ![Alt text](pics/table.png?raw=true "Optional Title")  
 
 * Output MS/MS spectrum with user-defined resolution (ppi) and Table 
 
## LICENSE  
### MIT
