import pandas as pd

def parse_files(cfg):
    inf = parse_input(cfg.input_file)
    outf = parse_output(cfg.output_file)
    return [inf, outf]

def parse_input(file):
    colums=['seqid', 'length', 'taxid', 'specialization']
    types={'seqid': 'str', 'length': 'uint64', 'taxid': 'str', 'specialization': 'str'}
    return pd.read_table(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_output(file):
    colums=['seqid', 'seqstart', 'seqend', 'length', 'taxid', 'binid', 'specialization']
    types={'seqid': 'str', 'start': 'uint64', 'end': 'uint64', 'length': 'uint64', 'taxid': 'str', 'binid': 'uint64', 'specialization': 'str'}
    return pd.read_table(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def sanity_check(cfg, inf, outf, missing_entries: int=0):
    if inf.empty: 
        print("Input file is empty")
        return False
    if outf.empty: 
        print("Output file is empty")
        return False

    if not missing_entries:
        if cfg.overlap_len:
            # Evaluate seqid and lengths separetely
            if not inf["seqid"].isin(outf["seqid"]).all():
                print("seqids not consistent between input and output")
                return False
            extra_length = (outf.shape[0] - inf.shape[0]) * cfg.overlap_len
            if (outf.length.sum()-extra_length != inf.length.sum()):
                print("lengths not consistent between input and output")
                return False
        else:
            # input grouped by seqid/length should match output
            if not inf.groupby(["seqid"]).sum()["length"].equals(outf.groupby(["seqid"]).sum()["length"]):
                print("seqids/lengths not consistent between input and output")
                return False
    else:
        # missing entries
        # Evaluate seqid and lengths separetely
        if inf["seqid"].unique().size != outf["seqid"].unique().size + missing_entries:
            print("seqids not consistent between input and output")
            return False
        if outf.length.sum() >= inf.length.sum():
            print("lengths not consistent between input and output")
            return False

    return True

class Config():

    def __init__(self,
        bin_exclusive: str=None, 
        bin_len: int=0, 
        bins: int=0, 
        fragment_len: int=0, 
        input_file: str=None,
        merged_file: str=None,
        nodes_file: str=None,
        overlap_len: int=0,
        output_file: str=None,
        pre_cluster: str=None,
        specialization: str=None,
        update_file: str=None,
        silent: bool=True,
        allow_merge: bool=False):

        self.bin_exclusive=bin_exclusive
        self.bin_len=bin_len
        self.bins=bins
        self.fragment_len=fragment_len
        self.input_file=input_file
        self.merged_file=merged_file
        self.nodes_file=nodes_file
        self.overlap_len=overlap_len
        self.output_file=output_file
        self.pre_cluster=pre_cluster
        self.specialization=specialization
        self.update_file=update_file
        self.silent=silent
        self.allow_merge=allow_merge