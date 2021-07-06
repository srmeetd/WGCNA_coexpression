from ruffus import *

import sys
import os
import sqlite3
import shutil
import CGATCore.Experiment as E
from CGATCore import Pipeline as P
import re
import glob
import collections
import CGAT.GTF as GTF
import CGATCore.IOTools as IOTools
import CGAT.BamTools.bamtools as BamTools
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC


# load options from the config file
P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS




# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions


# Specific pipeline tasks

@follows(mkdir("Coexpressed_modules.dir"))
@transform("dataset.dir/*", formatter("(.+).txt"), ["Coexpressed_modules.dir/{basename[0]}_all_modules.tsv", "Coexpressed_modules.dir/{basename[0]}_trib1_module.tsv"])
def WGCNA(infiles, outfiles):
#    job_memory = "10G"
    infile = infiles
    all_modules, trib1_module = outfiles
    output = P.snip(all_modules, "_all_modules.tsv")
    statement = '''Rscript /data/md1srd/Software/Co_expression_pipeline/run_wgcna.R --i %(infile)s --o %(output)s'''     
    
    P.run(statement)

@follows(WGCNA)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))




