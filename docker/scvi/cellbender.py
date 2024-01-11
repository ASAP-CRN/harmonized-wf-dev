import os
import argparse
import subprocess


parser = argparse.ArgumentParser(prog='cellbender-py')

parser.add_argument('-i', '--input', type=str)
parser.add_argument('-o', '--output', type=str)
parser.add_argument('-s', '--sample', type=str)
parser.add_argument('-f', '--fpr', type=str)

args = parser.parse_args()

os.chdir(os.path.join(os.getcwd(), 'data', args.sample))

cellbender_args = [
    '--input', args.input, 
    '--output', args.output,
    '--fpr', args.fpr,        
    '--cuda'
]


subprocess.run(['cellbender', 'remove-background'] + cellbender_args)


