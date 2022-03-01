# -*- coding: utf-8 -*-

"""
hmp2_workflows.tasks.file_conv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module contains functions that convert between different file types.

Copyright (c) 2017 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
"""

import os

import pandas as pd

from itertools import chain

from biobakery_workflows import utilities as bb_utils
from biobakery_workflows.tasks.sixteen_s import convert_to_biom_from_tsv


def deinterleave_fastq(workflow, input_files, output_dir, threads=1, compress=True):
    """Deinterleaves a FASTQ file producing paired-end FASTQ reads.

    Args:
        workflow (anadama2.Workflow): The AnADAMA2 Workflow object.
        input_files (list) A list of FASTQ files to deinterleave.
        output_dir (string): The output directory to write paired-end reads 
            too.
        compress (bool): Compress FASTQ files generated.
        threads (int): The number of threads/cores to be used if compressing
            paired ends reads.

    Requires:
        None

    Returns:
        list: A list of paired-end files.

    Example:
        from anadama2 import Workflow

        from hmp2_workflows.tasks.file_conv import deinterleave_fastq

        workflow = Workflow()
        paired_end_fastqs = deinterleave_fastq(workflow,
                                               ['foo.fasta', 'bar.fastq'],
                                               threads=4)

        print paired_end_fastqs                                                      
        # [[foo_R1.fastq.gz, foo_R2.fastq.gz], [bar_R1.fastq.gz, bar_R2.fastq.gz]]
    """
    paired_end_reads = []
    deinterleave_cmd = "deinterleave_fastq.sh < [depends[0]] [targets[0]] [targets[1]]"

    out_ext = "fastq"
    if compress:
        deinterleave_cmd += " " + "compress"
        out_ext = "fastq.gz"

    deinterleave_cmd += " " + str(threads)

    mate_1_files = bb_utils.name_files(map(os.path.basename, input_files),
                                        output_dir,
                                        tag="R1",
                                        extension=out_ext)
    mate_2_files = bb_utils.name_files(map(os.path.basename, input_files),
                                        output_dir,
                                        tag="R2",
                                        extension=out_ext)

    if "gz" in out_ext:
        mate_1_files = [fname.replace('.fastq_R1', '_R1.fastq') for fname in mate_1_files]
        mate_2_files = [fname.replace('.fastq_R2', '_R2.fastq') for fname in mate_2_files]
    output_files = zip(mate_1_files, mate_2_files)

    workflow.add_task_group_gridable(deinterleave_cmd,
                                     depends=input_files,
                                     targets=output_files,
                                     time=5*60,
                                     mem=4096,
                                     cores=threads)

    return output_files                                     


def bam_to_fastq(workflow, input_files, output_dir, paired_end=False,
                 compress=True, threads=1):
    """Converts BAM sequence files to a single interleaved FASTQ file using
    the samtools bam2fq utility.

    Args:
        workflow (anadama2.Workflow): The AnADAMA2 Workflow object to append 
            the BAM to FASTQ conversion step to.
        input_files (list): A list containing all BAM files to be converted.
        output_dir (string): The output directory to write converted files too.
        paired_end (bool): If True generated paired end files.
        compress (bool): Compress fastq files generated by samtools. 
        threads (int): The number of threads/cores to use for BAM -> FASTQ 
            conversion.

    Requires:
        bedtools 2.17+

    Returns:
        list: A list of the newly-converted FASTQ files.

    Example:
        from anadama2 import Workflow

        from hmp2_workflows.tasks.file_conv import bam_to_fastq

 
        workflow = Workflow()
        fastq_files = bam_to_fastq(workflow,
                                   ['/tmp/fooA.bam', '/tmp/fooB.bam'],
                                   '/seq/ibdmdbd/out_dir'])
    """
    sample_names = bb_utils.sample_names(input_files, '.bam')
    sorted_bams = bb_utils.name_files(sample_names, 
                                      output_dir, 
                                      subfolder="sort",
                                      tag="sorted", 
                                      extension="bam",
                                      create_folder=True)

    ## Gotta make sure our BAM file is sorted first
    workflow.add_task_group_gridable('sambamba sort -n -t [args[0]] -m 4GB -o [targets[0]] [depends[0]]',
                                     depends=input_files,
                                     targets=[os.path.splitext(bam)[0] for bam in sorted_bams],
                                     args=[threads],
                                     time=30*60,
                                     cores=threads,
                                     mem=4098)
    
    reformat_cmd = ("reformat.sh t=[args[0]] in=[depends[0]] out=stdout.fq primaryonly | " 
                    "reformat.sh t=[args[0]] in=stdin.fq out1=[targets[0]] ")
    if paired_end:
        mate_1_files = bb_utils.name_files(map(os.path.basename, input_files),
                                           output_dir,
                                           tag="R1",
                                           subfolder="fastq",
                                           extension="fastq",
                                           create_folder=True)
        mate_2_files = bb_utils.name_files(map(os.path.basename, input_files),
                                           output_dir,
                                           tag="R2",
                                           subfolder="fastq",
                                           extension="fastq",
                                           create_folder=True)

        mate_1_files = [fname.replace('.fastq_R1', '_R1.fastq') for fname in mate_1_files]
        mate_2_files = [fname.replace('.fastq_R2', '_R2.fastq') for fname in mate_2_files]
        output_files = zip(mate_1_files, mate_2_files)
        reformat_cmd += "out2=[targets[1]] "
    else:
        output_files = bb_utils.name_files(map(os.path.basename, input_files),
                                           output_dir,
                                           extension=".fastq")

    reformat_cmd += "interleaved addslash=t spaceslash=f"
    workflow.add_task_group_gridable(reformat_cmd,
                                     depends=input_files,
                                     targets=output_files,
                                     args=[threads],
                                     cores=threads,
                                     time=20*60,
                                     mem=4098)

    fastq_files = list(chain.from_iterable(output_files)) if paired_end else output_files

    if compress:
        fastq_files_compress = ["%s.gz" % fastq_file for fastq_file in fastq_files]

        workflow.add_task_group_gridable("pigz --best -p [args[0]] [depends[0]]",
                                         depends=fastq_files,
                                         targets=fastq_files_compress,
                                         args=[threads],
                                         cores=threads,
                                         time=10*60,
                                         mem=4098)
        fastq_files = fastq_files_compress
    
        workflow.add_task_group("rm -rf [targets[0]]",
                                targets=sorted_bams,
                                depends=fastq_files_compress)

    return fastq_files


def excel_to_csv(workflow, input_files, output_dir):
    """Converts an Excel file to a CSV file. Only attempts to convert the 
    first worksheet in the file and ignores the rest.

    Args:
        workflow (anadama2.Workflow): The AnADAMA2 workflow object.
        input_files (list): A list containing all Excel files to be converted.
        output_dir (string): The output directory to write converted CSV files
            too.

    Requires:
        None

    Returns:
        list: A list of newly-converted CSV files.
    """
    output_files = bb_utils.name_files(map(os.path.basename, input_files),
                                       output_dir,
                                       extension='csv')

    def _convert_excel_csv(task):
        """Helper function passed to AnADAMA2 doing the lifting of converting 
        the supplied Excel file to a CSV file using the pandas python library.
        """                                      
        excel_file = task.depends[0].name
        csv_out_file = task.targets[0].name

        excel_df = pd.read_excel(excel_file)
        excel_df.to_csv(csv_out_file)

    workflow.add_task_group(_convert_excel_csv,
                            depends=input_files,
                            targets=output_files)                                       

    return output_files                            


def batch_convert_tsv_to_biom(workflow, tsv_files): 
    """Batch converts tsv files to the biom format. BIOM files will be 
    deposited in the same folder as source TSV files and will carry the 
    same filenames.

    Args:
        workflow (anadama2.Workflow): The workflow object.
        tsv_files (list): A list containing all TSV files to be converted 
            to BIOM format.
    
    Requires:
        Biom v2: A tool for general use formatting of biological data.

    Returns: 
        list: A list containing paths to all converted BIOM files.

    Example:
        from anadama2 import Workflow
        from hmp2_workflows.tasks import common

        workflow = anadama2.Workflow()

        tsv_files = ['/tmp/foo.tsv', '/tmp/bar.tsv', '/tmp/baz.tsv']
        biom_files = common.batch_convert_tsv_to_biom(workflow, tsv_files)

        print biom_files
        ## ['/tmp/foo.biom', '/tmp/bar.biom', '/tmp/baz.biom']
    """
    biom_files = []

    tsv_fnames = bb_utils.sample_names(tsv_files, '.tsv')
    tsv_dir = os.path.dirname(tsv_files[0])

    biom_dir = os.path.join(tsv_dir, 'biom')
    bb_utils.create_folders(biom_dir)

    biom_files = [os.path.join(biom_dir, biom_fname) for biom_fname in 
                  bb_utils.name_files(tsv_fnames, biom_dir, extension='biom')]

    for (tsv_file, biom_file) in zip(tsv_files, biom_files):
        convert_to_biom_from_tsv(workflow, tsv_file, biom_file)

    return biom_files


def fix_CMMR_OTU_table_taxonomy_labels(workflow, otu_table, output_dir):
    """Takes an OTU table generated by CMMR and formats the table to make 
    sure it can be processed by some of the biobakery 16S viz functions.

    Args:
        workflow (anadama2.Workflow): The workflow object.
        otu_table (string): Path the CMMR OTU table
        output_dir (string): Directory

    Requires:
        None

    Returns:
        string: The path to the modified OTU table.
    """
    otu_basename = os.path.splitext(os.path.basename(otu_table))[0]
    fixed_otu_table = os.path.join(output_dir, otu_basename + "_taxonomy_fix.tsv")
    tax_level_label = {
        0: 'k',
        1: 'p',
        2: 'c',
        3: 'o',
        4: 'f',
        5: 'g'
    }

    def _label_taxonomy_levels(task):
        """Quick function that attempts to label the taxonomic functions in our OTU
        table in the format that some of the biobakery viz functions can handle.
        """
        fixed_otu_fh = open(fixed_otu_table, 'w')

        with open(otu_table) as otu_fh:
            for line in otu_fh:
                if line.startswith('#'):
                    fixed_otu_fh.write(line)
                    continue

                otu_elts = line.strip().split('\t')
                taxonomy = otu_elts[-1]
                new_taxonomy = []

                for (idx, tax_level) in enumerate(taxonomy.split('; ')):
                    tax_label = tax_level_label[idx]

                    if len(tax_level) == 3:
                        level_labeled = tax_level[-1::-1]
                    else:
                        if "__" in tax_level:
                            level_labeled = tax_label + tax_level
                        else:
                            level_labeled = "%s__%s" % (tax_label, tax_level)

                    new_taxonomy.append(level_labeled)

                fixed_otu_fh.write('\t'.join(otu_elts[:-1]) + '\t' + '; '.join(new_taxonomy) + '\n')

        fixed_otu_fh.close()
                                
    workflow.add_task(_label_taxonomy_levels,
                      depends=otu_table,
                      target= fixed_otu_table)

    return fixed_otu_table