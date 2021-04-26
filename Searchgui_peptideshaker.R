library("PGA")
Specta_files_directory="/media/DSRG4new/Hari/Spectra"
f<-list.files(pattern="*Trinity.fasta$")
setwd(Specta_files_directory)
system(find ./ -name '*.mzML' -exec msconvert --filter "peakPicking true 2-" --mgf {} \;)
for (k in 1:length(f))
{
i=f[k]
outdb <- createProDB4DenovoRNASeq(infa=i,outfile_name =i)
system(paste("awk 'BEGIN{ RS = \">\"; } { if ($0 !~ /#REV#/) { printf \">\"$0; } }' ",i,"_txFinder.fasta > ",i,"_txFinder_rev_removed1.fasta",sep='')) 
system(paste("awk 'BEGIN{FS=\"|\"}{if(/^>/){print \">\"$2}else{print $0}}' ",i,"_txFinder_rev_removed1.fasta >  ",i,"_txFinder_rev_removed_fasta_trimmed1.fasta",sep=''))
system(paste("cat reference_refseq.fasta ",i,"_txFinder_rev_removed_fasta_trimmed.fasta > ",i,"_exp_fasta_for_searching.fasta",sep=''))
system(paste("searchgui eu.isas.searchgui.cmd.FastaCLI -in ",i,"_exp_fasta_for_searching.fasta -decoy",sep=''))
system(paste("searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI -out ",i,".par -db   ",i,"_exp_fasta_for_searching.concatenated_target_decoy.fasta -frag_tol 0.05 -fixed_mods \"iTRAQ 4-plex of K,iTRAQ 4-plex of peptide N-term,Carbamidomethylation of C\" -variable_mods \"Acetylation of protein N-term,Deamidation of N,Oxidation of M\" -msgf_num_matches 1",sep=''))
system(paste("searchgui eu.isas.searchgui.cmd.SearchCLI -spectrum_files ",Specta_files_directory," -id_params ",i,".par -output_folder ", Specta_files_directory," -xtandem 1 -msgf 1 -tide 1  -output_default_name ",i,"_searchgui.out",sep=''))
}
f<-list.files(pattern="*zip")
f1<-list.files(pattern="*mgf$")
for (k in 1:length(f))
{
i<-f[k]
j<-f1[k]
system(paste("peptide-shaker eu.isas.peptideshaker.cmd.PeptideShakerCLI -experiment ",i," -sample ",i," -replicate 1 -identification_files ",Specta_files_directory,"/",i," -spectrum_files ",Specta_files_directory,"/",j," -out ",Specta_files_directory,"/",i,".cpsx",sep=''))
}
f<-list.files(pattern="*.cpsx$")
for(k in 1:length(f))
{
i=f[k]
system(paste("peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in ",Specta_files_directory,"/",i," -out_reports ",Specta_files_directory," -reports 3,6,9",sep=''))
}
