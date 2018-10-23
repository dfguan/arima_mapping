#!/usr/bin/env python3
import sys, os, json
sys.path.append("./runner")
from hpc import hpc
from manager import manager 
from multiprocessing import Process, Pool

def checkf(p):
    return os.path.isfile(p)

def checkd(p):
    return os.path.isdir(p)

def mkdir(d):
    os.makedirs(d)

def getfn(p):
    return os.path.basename(p)
def getd(p):
    dirn = os.path.dirname(p)
    return dirn if dirn else "."

def map_and_flt(para):
    [man, flter, ref, read_fl, out_fn] = para
    fn = getfn(out_fn)
    jcmd = "bwa mem -t 12 -B 8 {0} {1} | samtools view -Sb - | perl {2} | samtools view -Sb - > {3}".format(ref, read_fl, flter, out_fn)
    jout = "mapping_{}.o".format(fn)
    jerr = "mapping_{}.e".format(fn)
    
    p = hpc("lsf", cmd = jcmd, mem=20000, core=12, out=jout, err=jerr)
    return man.start([p])


def combine(para):
    [man, combiner, fn_lst, faidx_fn, maq, out_fn] = para
    fn = getfn(out_fn)
    jcmd = "perl {0} {1[0]} {1[1]} samtools {2} | samtools view -bS -t {3} - | samtools sort -o {4} -".format(combiner, fn_list, maq, faidx_fn, out_fn)
    jout = "combine_{}.o".format(fn)
    jerr = "combine_{}.e".format(fn)
    p = hpc("lsf", cmd = jcmd, mem=20000, core=12, out=jout, err=jerr)
    return man.start([p])
def add_rg(para):
    [man, picard, bam_fn, out_fn, rgid, rglb, rgsm] = para
    fn = getfn(out_fn)
    jcmd = "java -Xmx6g -jar {0} AddOrReplaceReadGroups INPUT={1} OUTPUT={2} ID={3} LB={4} SM={5} PL=ILLUMINA PU=none".format(picard, bam_fn, out_fn, rgid, rglb, rgsm)
    jout = "addrg_{}.o".format(fn)
    jerr = "addrg_{}.e".format(fn)
    p = hpc("lsf", cmd = jcmd, mem=8000, core=1, out=jout, err=jerr)
    return man.start([p])

def mark_dup(para):
    [picard, bam_fn, out_fn, matrix_fn, tmp_dir] = para
    fn = getfn(out_fn)
    jcmd = "java -Xms24G -XX:-UseGCOverheadLimit -Xmx24G -jar {0} MarkDuplicates INPUT={1} OUTPUT={2} METRICS_FILE={3} TMP_DIR={4} ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE".format(picard, bam_fn, out_fn, matrix_fn, tmp_dir)
    jout = "mkdup_{}.o".format(fn)
    jerr = "mkdup_{}.e".format(fn)
    p = hpc("lsf", cmd = jcmd, mem=30000, core=1, out=jout, err=jerr)
    return man.start([p])

def idx_fa(para):
    [man, ref] = para
    fn = getfn(ref)
    jcmd = ["samtools", "faidx", ref]
    jout = "idx_fa_{}.o"
    jout = "idx_fa_{}.e"
    p = hpc("bash", cmd = jcmd, out=jout, err=jerr)
    return man.start([p])

def idx_bam(para):
    [man, bam_fn] = para 
    fn = getfn(bam_fn)
    jcmd = ["samtools", "index", bam_fn]
    jout = "idx_bam_{}.o"
    jout = "idx_bam_{}.e"
    p = hpc("bash", cmd = jcmd, out=jout, err=jerr)
    return man.start([p])
def bwa_idx_fa(para):
    [man, ref] = para
    fn = getfn(bam_fn)
    jcmd = "bwa index {}".format(ref)
    jout = "bwa_idx_{}.o"
    jout = "bwa_idx_{}.e"
    p = hpc("lsf", cmd = jcmd, mem=10000, out=jout, err=jerr)
    return man.start([p])


def proc_pe(para):
    #mapp_flt
    [man, pe_fnlst, flter, combiner, picard, tmp_dir, maq, rgid, rglb, rgsm, ref, out_dir ] = para
    params = []
    out_bams = []
    for i in range(len(pe_fnlst)):
        out_fn = "{0}/{1}.maf.bam".format(out_dir, getfn(pe_fnlst[0]))
        out_bams.append(out_fn)
        params.append([man, flter, ref, pe_fnlst[i], out_fn])
    pl = Pool(processes=len(pe_fnlst))
    rtvs = pl.map(map_and_flt, params)
    suc = all(v == 0 for v in rtvs)
    rtv = 0 if suc else 1
    if not rtv:
        pl = Pool(process=1)
        in_bams = out_bams
        out_bam = "{0}/{1}.cmb.bam".format(out_dir, get_rm_prefix(out_bams[0]))
        params = [[man, combiner, in_bams, "{}.fai".format(ref), maq, out_bam]]
        rtv = pl.map(combine, params)
        if not rtv:
            in_bam = out_bam
            out_bam = "{0}/{1}.rg.bam".format(out_dir, get_rm_prefix(in_bam))
            pl = Pool(process=1)
            params = [[man, picard, in_bam, out_bam,rgid, rglb, rgsm]]
            rtv = pl.map(add_rg, params)
            if not rtv:
                in_bam = out_bam
                out_bam = "{0}/{1}.mkdup.bam".format(out_dir, get_rm_prefix(in_bam))
                out_matrix = "{0}/{1}.mkdup.matrix.txt".format(out_dir, get_rm_prefix(in_bam))
                pl = Pool(process=1)

                params = [[man, picard, in_bam, out_bam, out_matrix, tmp_dir]]
                rtv = pl.map(mark_dup, params)
                if not rtv: 
                    in_bam = out_bam
                    pl = Pool(process=1)
                    params = [[man, in_bam]]
                    rtv = pl.map(idx_bam, params)
    return rtv
if __name__== "__main__":
    
    if len(sys.argv) < 2:
        print ("arima_mapping <config>")
        sys.exit(1)
    ma = manager("sys.config", retries=2) 
    config_fn = sys.argv[1]
    f = open(config_fn, "r")
    cfg_dict = json.load(f)
    fofn = cfg_dict["fofn"]
    hic_fn_lst = []

    ref = cfg_dict["REF"]
    fai_fn = "{}.fai".format(ref)
    if not checkf(fai_fn):
        pl = Pool(processes=1)
        rtv = pl.map(idx_fa, [[man, ref]])
        if rtv:
            print ("Failed to generate faidx for reference")
            sys.exit(1)
    bwa_idx = "{}.bwt".format(ref)
    if not checkf(bwa_idx):
        pl = Pool(processes=1)
        rtv = pl.map(bwa_idx_fa, [[man, ref]])
        if rtv:
            print ("Failed to generate bwa index for reference")
            sys.exit(1)
    with open(fofn, "r") as f:
        for ln in f:
            hic_fn_lst.append(ln.strip('\n').split('\t'))
        f.close()
    params = []
    for i in range(len(hic_fn_lst)):
        params.append([man, hic_fn_lst[i], cfg_dict["FILTER"], cfg_dict["COMBINER"], cfg_dict["PICARD"], cfg_dict["TMP_DIR"],cfg_dict["MAPQ_FILTER"], cfg_dict["SRA"], cfg_dict["SRA"], cfg_dict["LABEL"], ref, cfg_dict["OUT_DIR"]])
    pl = Pool(processes = len(hic_fn_lst))
    rtv = pl.map(proc_pe, params)
    if not rtv:
        print("command run successfully")




