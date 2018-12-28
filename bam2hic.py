#env python3
# Author: Dengfeng Guan
# Aim   : convert sam/bam to hic
# Date  : 2018/12/28
# Dependencies: pysam, runner

import os, sys, argparse
import pysam
from runner.runner.manager import manager
from runner.runner.hpc import hpc

def sf2hic(man, jboxt_path, in_fn, out_fn, chr_bed, out_dir):
    jcmd = "java -Xmx2g -jar {3} pre {0} {1} {2}".format(in_fn, out_fn, chr_bed, jboxt_path).split(' ')
    jout = "{}/jboxt.o".format(out_dir)
    jerr = "{}/jobxt.e".format(out_dir) 
    # j = hpc("lfs", cmd = jcmd, out=jout, err=jerr, core = 1, queue=long, mem = 20000)
    j = hpc(cmd = jcmd, out=jout, err=jerr, core = 1, queue="long", mem = 20000)
    return man.start([j])

#sort according to chromosome names
def sort_by_chrn(in_fn, out_fn, out_dir):
    jcmd = "sort -k2,2 -k6,6 {0}".format(in_fn).split(' ')
    jout = out_fn 
    jerr = "{}/srt.e".format(out_dir) 
    # j = hpc("lfs", cmd = jcmd, out=jout, err=jerr, core = 1, queue=long, mem = 20000)
    j = hpc(cmd = jcmd, out=jout, err=jerr, core = 1, queue="long", mem = 20000)
    return man.start([j])
    



#
def process_bams(bam_fns, out_fn, chr_bed):
    hasw_hdr = False
    with open(out_fn, "w") as fout:
        for bam_fn in bam_fns:
            bf = pysam.AlignmentFile(bam_fn, "rb")
            if not hasw_hdr:
                with open(chr_bed, "w") as fchr:
                    for ele in bf.header["SQ"]:
                        fchr.write("{0}\t{1}\n".format(ele['SN'], ele['LN']))
                    fchr.close()
                hasw_hdr = True
            qn = ""
            fv = []
            for read in bf:
                cur_qn = read.query_name
                if cur_qn == qn:
                    if not read.is_unmapped:
                        c1 = read.cigartuples[0]
                        cl = read.cigartuples[-1]
                        if (read.is_reverse and cl[0] == 0) or (not read.is_reverse and c1[0] == 0):
                            fv.append(read)
                else:
                    if qn != "":
                        if len(fv) == 2 and fv[0].is_read1 != fv[1].is_read1 and fv[0].mapping_quality > 10 and fv[1].mapping_quality > 10: #not the same read 
                            if fv[0].reference_name > fv[1].reference_name:
                                small = 1
                            else:
                                small = 0
                            r1_details = [1 if fv[small].is_reverse else 0, fv[small].reference_name, fv[small].reference_start, 8]
                            r2_details = [1 if fv[1-small].is_reverse else 0, fv[1-small].reference_name, fv[1-small].reference_start, 16]
                            
                            fout.write("{0[0]}\t{0[1]}\t{0[2]}\t{0[3]}\t{1[0]}\t{1[1]}\t{1[2]}\t{1[3]}\n".format(r1_details, r2_details))
                        
                    qn = cur_qn
                    if not read.is_unmapped:
                        c1 = read.cigartuples[0]
                        cl = read.cigartuples[-1]
                        if (read.is_reverse and cl[0] == 0) or (not read.is_reverse and c1[0] == 0):
                            fv = [read]
        fout.close() 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='produce hic file for juicebox with bam files')
    parser.add_argument('-p', help='Prefix of output files',dest='prefix', default="st")
    parser.add_argument('-b', help='Path of juicebox tool', dest='jboxt_path', default="./")
    parser.add_argument('bams', help="bam files", nargs='+')
    parser.add_argument('out_dir', help="output directory")

    args = parser.parse_args()    
    
    man = manager(retries=2, wait=2)
    
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    
    out_tmp = "{0}/{1}.sft.bed".format(args.out_dir, args.prefix)
    out_chr_bed = "{0}/{1}.chrs.sz.bed".format(args.out_dir, args.prefix)
    process_bams(args.bams, out_tmp, out_chr_bed)
    print ("finished generating short format contact map with bam files") 

    in_fn = out_tmp 
    out_srt_tmp = "{0}/{1}.sft.srt.bed".format(args.out_dir, args.prefix)
    rtn = sort_by_chrn(in_fn, out_srt_tmp, args.out_dir)
    if not rtn: 
        print ("finished sorting contact map") 
    else:
        print ("failed to sort contact map") 
        sys.exit(1)
    in_fn = out_srt_tmp
    in_chr_bed = out_chr_bed
    out_fn = "{0}/{1}.hic".format(args.out_dir, args.prefix)
    rtn = sf2hic(man, args.jboxt_path, in_fn, out_fn, out_chr_bed, args.out_dir)
    
    if not rtn: 
        print ("finished generating hic file") 
    else:
        print ("failed to generate hic file") 
    
    
