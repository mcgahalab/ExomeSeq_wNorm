import glob
import pandas as pd
from snakemake.utils import validate
#from snakemake.remote import FTP
#ftp = FTP.RemoteProvider()
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(
    config["sample_file"]
).set_index(
    "sample", drop=False
)

def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired

def get_map_reads_input_R1(wildcards):
    if not is_activated("mergeReads"):
        if config["trimming"]["activate"]:
            return expand(
                "results/trimmed/{sample}_{unit}_R1.fastq.gz",
                unit=units.loc[wildcards.sample, "unit_name"],
                sample=wildcards.sample,
            )
        unit = units.loc[wildcards.sample]
        if all(pd.isna(unit["fq1"])):
            # SRA sample (always paired-end for now)
            accession = unit["sra"]
            return expand("sra/{accession}_R1.fastq", accession=accession)
        sample_units = units.loc[wildcards.sample]
        return sample_units["fq1"]
    if is_paired_end(wildcards.sample):
        return "results/merged/{sample}_R1.fastq.gz"
    return "results/merged/{sample}_single.fastq.gz"


def get_map_reads_input_R2(wildcards):
    if is_paired_end(wildcards.sample):
        if not is_activated("mergeReads"):
            if config["trimming"]["activate"]:
                return expand(
                    "results/trimmed/{sample}_{unit}_R1.fastq.gz",
                    unit=units.loc[wildcards.sample, "unit_name"],
                    sample=wildcards.sample,
                )
            unit = units.loc[wildcards.sample]
            if all(pd.isna(unit["fq1"])):
                # SRA sample (always paired-end for now)
                accession = unit["sra"]
                return expand("sra/{accession}_R2.fastq", accession=accession)
            sample_units = units.loc[wildcards.sample]
            return sample_units["fq2"]
        return ("results/merged/{sample}_R2.fastq.gz",)
    return ""


def get_r1(wildcards):
    return samples.read1[wildcards.sample]

def get_r2(wildcards):
    return samples.read2[wildcards.sample]

def is_control(wildcards):
    control = samples.loc[ (wildcards.sample), "isControl" ]
    return pd.isna(control) or pd.isnull(control)

def has_a_control(wildcards):
    # Assumes samples are in [SAMPLE]_[CELLTYPE]_[ANTIBODY] format
    control = samples.loc[ (wildcards.sample), "control" ]
    if pd.isna(control) or pd.isnull(control):
        return False
    else:
        return control != "_".join([wildcards.sample])

def get_sample_control(wildcards, retself=True):
    # Assumes control is a [SAMPLE]_[CELLTYPE]_[ANTIBODY] format
    control = samples.loc[ (wildcards.sample), "control" ]
    if pd.isna(control) or pd.isnull(control):
        control = [ wildcards.sample ] if retself else ""
    else:
        control = control.split("_")
    return "_".join(control)

intervals = pd.read_table(
    config["bed_file"]
).set_index(
    "interval", drop=False
)
def get_intervals(wildcards):
    inter = wildcards.interval
    bed = "/cluster/projects/cesconlab/Projects/GENIUS/TGL10_legacy_data/ExomeSeqNorm/resources/hg38_bedV5/" + inter + ".bed"
    return bed

def get_MuTect2_output(wildcards):
    res = []
    for i in intervals.itertuples():
        res.append(
            "results/MuTect2/{}/{}_{}.mut2.vcf".format(
                wildcards.sample, wildcards.sample, i.interval
            )
        )
    return res

indel_vcfs = pd.read_table(
    config["indel_vcf"]
,dtype={'indel': object}).set_index(
    "indel", drop=False
)
def get_indels(wildcards):
    inter = wildcards.indel
    indel = str(inter) + "_hg38" + ".vcf"
    #indel = str(inter) + ".vcf"
    return indel

snv_vcfs = pd.read_table(
    config["snv_vcf"]
,dtype={'snv': object}).set_index(
    "snv", drop=False
)

def get_snvs(wildcards):
    inter = wildcards.snv
    #snv = str(inter) + "_hg38" + ".vcf"
    snv = str(inter) + "_hg38" + ".vcf"
    return snv

def get_maf_output(wildcards, type='indel'):
    res = []
    if type == 'snv':
        for v in snv_vcfs.itertuples():
            res.append(
                "results/MAF_38_final/{}/{}/{}.maf".format(
                    type, wildcards.sample, v.snv
                )
            )
    elif type == 'indel':
        for v in indel_vcfs.itertuples():
            res.append(
                "results/MAF_38_final/{}/{}/{}.maf".format(
                    type, wildcards.sample, v.indel
                )
            )
    return res

def get_snv_intersects(wildcards):
    return config["intersects"]["snv"][wildcards.snv]

def get_indel_intersects(wildcards):
    return config["intersects"]["indel"][wildcards.snv]
