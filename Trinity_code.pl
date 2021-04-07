#perl Trinity.pl Fastq_files_directory

#perl code 
use strict;
use diagnostics;

#pwd should be the one where the files are downloaded; Open command prompt in that directory and use the following command to run the script
#perl RNASeq.pl /media/DSRG4new/3D or any other path in which your files are stored

my $dir = $ARGV[0];
my $file;
my $file1;
my $file2;
my @readFiles;
my $count;
my $outFile;


opendir(DIR, $dir) or die $!;
    while ($file = readdir(DIR)) {
        # Use a regular expression to ignore files beginning with a period
        #next if ($file =~ m/(^\.|_R2.fastq$|out$|bam$|gz$|gz.1$|pl$|sh$)/);
	next if ( -d $file );
	if($file =~ /(.*)\_R1\.fastq$/){
		print "$file\n";		
		@readFiles = split(/\_R1/,$file);
		$file1 = $dir."/".$readFiles[0]."_R1.fastq";
		print "$file1\n";
		$file2 = $dir."/".$readFiles[0]."_R2.fastq";
		print "$file2\n";
		$outFile = $readFiles[0]."OUT";
		print "$readFiles[0]\n";
                `Trinity --seqType fq --normalize_by_read_set --left $file1.fastq --right $file2.fastq --trimmomatic --full_cleanup --CPU 30 --max_memory 50G --bflyCPU 10 --bflyHeapSpaceMax 4G --output Trinity.$file --monitoring --verbose`;		
		$count++; 
	}
    }	
#
