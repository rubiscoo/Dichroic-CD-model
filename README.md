# Ensemble dichroic model peptide helicity estimator

Python script to estimate helical peptide helicity from measured molar residue ellipticity (MRE).

<b> Requirements: </b>
- installed python3
- standard libraries: numpy, scipy

<b> Instructions for use: </b>
1. Clone/download the repository. The folder should contain the python script (.py) along with folders <i>Q_total</i> and <i>Q_double_H</i>, that contain partition function polynomials needed for calculation.
2. Navigate to the folder using terminal/command line [ cd ./folder/folder/ ]
3. Run the script using python3 command [ python3 script.py ]
4. Define three requested input parameters: <b>N</b> (number of residues), <b>temperature</b> (in °C) and measured <b>MRE</b> (in deg cm<sup>2</sup> dmol<sup>-1</sup> per peptide unit)
   
<b> OUTPUT (result): </b> Program returns propagation parameter <b><i>w</i></b> and fractional helicity <b><i>f<sub>H</sub></i></b>.

Default parameters are those obtained from global fit of CD data (AAKAA<sub>n</sub> series of peptides and CCs). See article "Estimation of peptide helicity from circular dichroism using the ensemble model ".
