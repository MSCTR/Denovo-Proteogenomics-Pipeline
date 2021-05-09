#!/usr/bin/python
import os, sys, getopt, re, subprocess

#python trinity_code.py /path/to/fastq/dir
def usage():
    print("""
trinity_code.py :   Passes paired end FASTQ files to Trinity for de novo assembly of RNA transcripts. 

trinity_code.py [-h] [-s <single end>] </path/to/directory>

    -h              Print this message

    -s              If single end reads, pass only one file to Trinity

    <directory>     Directory containg reads in FASTQ format
""")    

# Define shell arguments
opt, req = getopt.getopt(sys.argv[1:], 'hs')
opts = {}
for k,v in opt:
    opts[k] = v
if '-h' in opts.keys():
    usage(); sys.exit()
if '-s' in opts.keys():
    pass

dir_path = sys.argv[1]

# Ensure non-empty directory
if not os.listdir(dir_path):
    sys.exit("Directory %s is empty." %dir_path)

contents = os.listdir(dir_path)
r1_pattern = re.compile('(.*)_R1.fastq$') # search for fastq files ending in '_R1'
for i in range(len(contents)):
    if r1_pattern.match(contents[i]):
        id = contents[i].split('_R1')
        if id[0]+'_R2.fastq' in contents: # ensures only paired reads
            print(id[0])
            file1 = dir_path + "/" + id[0] + "_R1.fastq"
            print(file1)
            file2 = dir_path + "/" + id[0] + "_R2.fastq"
            print(file2)
            os.system(f"Trinity --seqType fq --normalize_by_read_set --left {file1} --right {file2} --trimmomatic --full_cleanup --CPU 30 --max_memory 50G --bflyCPU 10 --bflyHeapSpaceMax 4G --output {dir_path}/Trinity.{id[0]} --monitoring --verbose")
            

    