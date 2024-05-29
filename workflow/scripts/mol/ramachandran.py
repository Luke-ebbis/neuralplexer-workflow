#!/usr/bin/env python
# -*- coding: utf-8 -*-
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys 
import numpy as np
import itertools
import logging
import polars as pl
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral, Ramachandran
import matplotlib.pyplot as plt

# Configure the logging system
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filename='app.log',
                    filemode='w')

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(
    logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))

# Plotting Ramachandran analysis using MD analysis
# https://docs.mdanalysis.org/1.1.0/documentation_pages/analysis/dihedrals.html
# OR ---
# https://mdtraj.org/1.9.4/examples/ramachandran-plot.html

def load_pdb(path: str):
    u = mda.Universe(path)
    return u

def load_pdb2(path: str):
    """Load a PDB structure

    Parameters
    ==========
    path:
        The path leading to the structure.

    Returns
    ======
    An md.Trajectory.
    """
    trajectory = md.load(path)
    return trajectory


def compute_ramachandran(universe):
    residues = universe.select_atoms(f"resid 1-{len(universe.residues)}")
    ram = Ramachandran(residues).run()
    return ram


def plot_ramachandran(ram):
    fig, ax = plt.subplots(figsize=plt.figaspect(1))
    ram.plot(ax)
    return fig, ax

def remove_first_last(df):
    return df.slice(1, df.height -2)

def ramachandran_to_df(ram):
    res = ram.atomgroup
    resids, names, segids = (pl.from_numpy(res.resids, schema=["resids"]),
                     pl.from_numpy(res.resnames, schema=["resnames"]), 
                     pl.from_numpy(res.segids, schema=["segids"])
                     )
    indexes = pl.concat([resids, names, segids],
                        how="horizontal").unique().sort("resids")
    meta = remove_first_last(indexes)
    print(meta)
    frame_df = []
    for frame_id, frame in enumerate(ram.results['angles']):
        analysis = pl.from_numpy(frame, schema=["psi", "phi"])

        current_frame = (pl.concat([meta, analysis],
                                  how="horizontal").
                         with_columns(pl.lit(frame_id).alias("frame_id"))
                         )
        print(current_frame)
        frame_df.append(current_frame)
    df = pl.concat(frame_df)
    return df


def write_plot(ramachandran, path):
    fig, ax = plot_ramachandran(ramachandran)
    fig.savefig(path)



def main():
    # Get the root logger
    input_structure = snakemake.input[0]
    output_plot = snakemake.output[0]
    output_analysis = snakemake.output[1]
    logger = logging.getLogger()
    logger.addHandler(console_handler)
    path = input_structure
    pdb = load_pdb(path)
    ramachandran = compute_ramachandran(pdb)

    # as a dataframe
    df = ramachandran_to_df(ramachandran)
    df.write_csv(output_analysis)
    
    # as a plot
    write_plot(ramachandran, output_plot)



if __name__ == "__main__":

    main()
