function [MasterSheet, Genes, Headers, GeneTable] = Import_Seq(directory,f_type,d_cols)
%  *                                                                     *
%  *              Inputs: directory = data directory                     *
%  *                      f_type = file type to search for in directory  *
%  *                      d_cols = [Genes Column, Count Column]          *
%  *              (Gene Column is the Column # only including text data) *
%  *              (Count Column is the Column # only including # data)   *
%  *                                                                     *
%  *              Outputs: MasterSheet = All count data                  *
%  *                       Genes = Vector of Gene IDs                    *
%  *                       Headers = Vector of Subject IDs               *
%  *                                                                     *
% Setting up Master Gene Table
cd(directory); % Change to data directory
Headers={};
MasterSheet=[];


% Import Gene Count data from an Excel Sheet
seq_files=dir(f_type); %find all gene tables using file extension
for j=1:length(seq_files) % loop through seq data files
    [tok]=strtok(seq_files(j).name, '.'); % Parses Animal Name From file Name
    disp(['Importing - ' tok]); % Display Importing File
    [a b]=strtok(tok); 
    b=b(~isspace(b)); % Find Spaces in File Name
    tok=[b a]; % Remove Spaces From file name
    Headers{j}=tok; % Add animal name to header array
    tmp=importdata(seq_files(j).name); % Import seq data
    
    counts=round(tmp.data(:,d_cols(2)));
    MasterSheet(:,j)=counts; % Gene counts (FKPM) for each Subject  ||IMPORTANT for "tmp.data(:,X)"  X=COLUMN WITH GENE COUNTS||
    if j==1
       Genes=tmp.textdata(2:end,d_cols(1)); % Get Gene Labels
       GeneTable=table(Genes,MasterSheet(:,j));
    else
       GeneTable=[GeneTable table(MasterSheet(:,j))];
       GeneTable.Properties.VariableNames=['Genes',Headers(1:j)];
    end
    
end
end

