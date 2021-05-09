library("PGA")
Spectra_files_directory="./"
f <- list.files(pattern="*Trinity.fasta$")
setwd(Spectra_files_directory)
#system(find ./ -name '*.mzML' -exec msconvert --filter "peakPicking true 2-" --mgf {} \;)

for (k in 1:length(f))
{
i=f[k]
outdb <- createProDB4DenovoRNASeq(infa=i,outfile_name =i)
system(paste("awk 'BEGIN{ RS = \">\"; } { if ($0 !~ /#REV#/) { printf \">\"$0; } }' ",i,"_txFinder.fasta > ",i,"_txFinder_rev_removed1.fasta",sep=''))
system(paste("awk 'BEGIN{FS=\"|\"}{if(/^>/){print \">\"$2}else{print $0}}' ",i,"_txFinder_rev_removed1.fasta >  ",i,"_txFinder_rev_removed_fasta_trimmed1.fasta",sep=''))
system(paste("cat ./resources/reference_hg19.fasta ",i,"_txFinder_rev_removed_fasta_trimmed1.fasta > ",i,"_exp_fasta_for_searching.fasta",sep=''))
system(paste("java -cp ./resources/searchgui/SearchGUI-3.3.17.jar eu.isas.searchgui.cmd.FastaCLI -in ",i,"_exp_fasta_for_searching.fasta -decoy",sep=''))
system(paste("sudo java -cp ./resources/searchgui/SearchGUI-3.3.17.jar eu.isas.searchgui.cmd.IdentificationParametersCLI -out ",i,".par -db   ",i,"_exp_fasta_for_searching_concatenated_target_decoy.fasta -frag_tol 0.05 -fixed_mods \"iTRAQ 8-plex of K,iTRAQ 8-plex of peptide N-term,Carbamidomethylation of C\" -variable_mods \"Acetylation of protein N-term,Deamidation of N,Oxidation of M\" -msgf_num_matches 1",sep=''))
system(paste("sudo java -cp ./resources/searchgui/SearchGUI-3.3.17.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files ",Spectra_files_directory," -id_params ",i,".par -output_folder ", Spectra_files_directory," -xtandem 1 -msgf 1 -tide 1  -output_default_name ",i,"_searchgui.out",sep=''))
}

f<-list.files(pattern="*zip")
f1<-list.files(pattern="*mgf$")
for (k in 1:length(f))
{
i<-f[k]
j<-f1[k]
system(paste("sudo java -cp ./resources/peptideshaker/PeptideShaker-1.16.42.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -experiment ",i," -sample ",i," -replicate 1 -identification_files ",Spectra_files_directory,"/",i," -spectrum_files ",Spectra_files_directory,"/",j," -out ",Spectra_files_directory,"/",i,".cpsx",sep=''))
}
f<-list.files(pattern="*.cpsx$")
for(k in 1:length(f))
{
i=f[k]
system(paste("java -cp ./resources/peptideshaker/PeptideShaker-1.16.42.jar eu.isas.peptideshaker.cmd.ReportCLI -in ",Spectra_files_directory,"/",i," -out_reports ",Spectra_files_directory," -reports 3,6,9",sep=''))
}
