
## Imports ##

from optparse import OptionParser
import datetime
import os
import math

## Import arguments ## 

parser = OptionParser()
parser.add_option("-I", "--Input_file", action="store", type="string", dest="Input_filename",
                  help="Input .BAdabouM file", metavar="File.BAdabouM")
parser.add_option("-r", "--Reference", action="store", type="string", dest="Ref",
                  help="Reference Used to align reads", metavar="Reference", default="Null")
parser.add_option("-O", "--Output_file", action="store", type="string", dest="Output_filename",
                  help="output .vcf file", metavar="File.vcf")
(options, args) = parser.parse_args()

## Transform .BAdabouM in .vcf ##

now = datetime.datetime.now()

## Open files 
Input = open(options.Input_filename, "r")
Output = open(options.Output_filename, "w")

## Write header

# File Format
Output.write("##fileformat=VCFv4.1\n")

# date
Output.write("##filedate=")
Output.write(now.strftime("%Y%m%d\n"))

# Reference
if options.Ref != "Null":
    Output.write("##reference=")
    Output.write(options.Ref)
    Output.write("\n")

# infos 
Output.write('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">\n')
Output.write('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">\n')
Output.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
Output.write('##INFO=<ID=SVLEN,Number=-1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
Output.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
Output.write('##ALT=<ID=DEL,Description="Deletion">\n')
Output.write('##ALT=<ID=INS,Description="Insertion of novel sequence">\n')
Output.write('##ALT=<ID=INV,Description="Inversion">\n')
Output.write('##ALT=<ID=CNV,Description="Copy number variable region">\n')
Output.write('##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">\n')

# colNames
Name=os.path.basename(str(options.Input_filename)).split('.')[0]
Output.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
Output.write(Name)
Output.write('\n')


## Process SVs

i=0
for line in Input.readlines():
    if (i==0): i=1
    else :
        Line = str(line.split('\n')[0])
        Line = line.split('\t')
        
        CHROM=Line[0]
        
        POS=math.floor((int(Line[1])+int(Line[2]))/2)

        ID="."
        REF="A"
        #ALT="N"
        QUAL="."
        FILTER="PASS"

        # INFO
        SVTYPE=Line[6]
        END=math.floor((int(Line[4])+int(Line[5]))/2)
        
        SVLEN=END-POS
        
        CIPOS_1=int(Line[1])-POS
        CIPOS_2=int(Line[2])-POS
        CIEND_1=int(Line[4])-END
        CIEND_2=int(Line[5])-END
        
        FORMAT="GT"
        
        GENO="1/1"
        
        Output.write(CHROM)
        Output.write("\t")
        
        Output.write(str(POS))
        Output.write("\t")
        
        Output.write(ID)
        Output.write("\t")
        
        Output.write(REF)
        Output.write("\t")
        
        Output.write("<")        
        Output.write(SVTYPE)
        Output.write(">\t")
        
        Output.write(QUAL)
        Output.write("\t")
        
        Output.write(FILTER)
        Output.write("\t")
        
        Output.write("IMPRECISE;SVTYPE=")
        Output.write(SVTYPE)
        Output.write(";END=")
        Output.write(str(END))
        Output.write(";SVLEN=")
        Output.write(str(SVLEN))
        Output.write(";CIPOS=")
        Output.write(str(CIPOS_1))
        Output.write(",")
        Output.write(str(CIPOS_2))
        Output.write(";CIEND=")
        Output.write(str(CIEND_1))
        Output.write(",")
        Output.write(str(CIEND_2))
        Output.write("\t")
        
        Output.write(FORMAT)
        Output.write("\t")
        
        Output.write(GENO)
        Output.write('\n')
