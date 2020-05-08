clear all

% Enter Directory of Salmon Gene Quantification for IP
directory='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Gene Quantification\IP';

% File Type and Data Columns [Gene TPM]
f_type='*.tabular';
d_cols = [1 3];

% Import Samples
[MasterSheet, Genes, Headers, GeneTable] = Import_Seq(directory,f_type,d_cols);

% Set Data Folder and Fraction Name
cd('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\Data Sheets');
writetable(GeneTable,'GeneTableIP.csv')

% Enter Directory of Salmon Gene Quantification for Input
directory='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Gene Quantification\Input';

% File Type and Data Columns [Gene TPM]
f_type='*.tabular';
d_cols = [1 3];

% Import Samples
[MasterSheet, Genes, Headers, GeneTable] = Import_Seq(directory,f_type,d_cols);

% Set Data Folder and Fraction Name
cd('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\Data Sheets');
writetable(GeneTable,'GeneTableIN.csv')

% Enter Directory of Salmon Gene Quantification for IP
directory='Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\Gene Quantification\IP_Minus';

% File Type and Data Columns [Gene TPM]
f_type='*.tabular';
d_cols = [1 3];

% Import Samples
[MasterSheet, Genes, Headers, GeneTable] = Import_Seq(directory,f_type,d_cols);

% Set Data Folder and Fraction Name
cd('Z:\Neumaier Lab\Morphine Grant\RNA Seq 2\WGCNA for Publication\Data Sheets');
writetable(GeneTable,'GeneTableIP_Minus.csv')