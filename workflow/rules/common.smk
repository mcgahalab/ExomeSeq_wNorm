samples = pd.read_table(
    config["sample_file"]
).set_index(
    "sample", drop=False
)

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
    bed = "/cluster/home/selghamr/workflows/ExomeSeq_wNorm/resources/hg38_bedV5/" + inter + ".bed"
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
