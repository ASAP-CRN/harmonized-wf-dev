import os
import argparse
import subprocess
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Pre-preprocess')
# Add arguments
parser.add_argument('--raw-counts', dest='raw_counts', type=str, 
                    help='Unfiltered feature-barcode matrices HDF5 output by cellranger')
parser.add_argument('--output-name', dest='output', type=str, 
                    help='Output file to save cellbender outputs to object to')
parser.add_argument('--fpr', dest='fpr', type=str, 
                    help='False positive rate', default='0.0')

args = parser.parse_args()


cellbender_args = [
    '--cuda'
    '--input', args.raw_counts, 
    '--output', args.output,
    '--fpr', args.fpr,        

]

# generates _out
subprocess.run(['cellbender', 'remove-background'] + cellbender_args)

"""
https://cellbender.readthedocs.io/en/latest/usage/index.html
This command will produce nine output files:
    {output}output_report.html: 
      HTML report including plots and commentary, along with any warnings or suggestions for improved parameter settings.
    {output}output.h5: 
        Full count matrix as an h5 file, with background RNA removed. This file contains all the original droplet barcodes.
    {output}output_filtered.h5: 
        Filtered count matrix as an h5 file, with background RNA removed. The word “filtered” means that this file contains 
        only the droplets which were determined to have a > 50% posterior probability of containing cells.
    {output}output_cell_barcodes.csv: 
        CSV file containing all the droplet barcodes which were determined to have a > 50% posterior probability of 
        containing cells. Barcodes are written in plain text. This information is also contained in each of the above 
        outputs, but is included as a separate output for convenient use in certain downstream applications.
    {output}output.pdf: 
        PDF file that provides a standard graphical summary of the inference procedure.
    {output}output.log: 
        Log file produced by the cellbender remove-background run.
    {output}output_metrics.csv: 
        Metrics describing the run, potentially to be used to flag problematic runs when using CellBender as part of a l
        arge-scale automated pipeline.
    {output}ckpt.tar.gz: 
        Checkpoint file which contains the trained model and the full posterior.
    {output}output_posterior.h5: 
        The full posterior probability of noise counts. This is not normally used downstream.

"""
