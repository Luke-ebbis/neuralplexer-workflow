#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Wrapper to the neuralPlexer programme.
"""

import argparse
import json
import subprocess
from typing import Dict, List
from os import environ
from pathlib import Path
DESCRIPTION = "Wrapper to the neuralPlexer programme"

def parse_json(file_path: str) -> Dict:
    with open(file_path, 'r') as file:
        data = json.load(file)
        return data


def check_environment():
    if environ.get('PIXI_ENVIRONMENT_NAME') is None:
        raise EnvironmentError("Run this programme from a pixi setting.")


def concat(items: List) -> str:
    return "|".join(items)


def receptor(job: Dict) -> str:
    sequence = []
    for seq in job:
        for number in range(seq['proteinChain']['count']):
            sequence.append(seq['proteinChain']['sequence'])
    seq_string = concat(sequence)
    return seq_string


def ligand(job: Dict) -> str:
    ligand = []
    for lig in job:
        for number in range(lig['molecule']['count']):
            ligand.append(lig['molecule']['sdf'])
    lig_string = concat(ligand)
    return lig_string


def make_jobs(job_data: Dict,
              image: str,
              checkpoint: str,
              output_path: str,
              cuda=True) -> List:
    """

    usage: neuralplexer-inference [-h] --task TASK [--sample-id SAMPLE_ID]
                                  [--template-id TEMPLATE_ID] [--cuda]
                                  [--model-checkpoint MODEL_CHECKPOINT]
                                  [--input-ligand INPUT_LIGAND]
                                  [--input-receptor INPUT_RECEPTOR]
                                  [--input-template INPUT_TEMPLATE]
                                  [--out-path OUT_PATH] [--n-samples N_SAMPLES]
                                  [--chunk-size CHUNK_SIZE]
                                  [--num-steps NUM_STEPS]
                                  [--latent-model LATENT_MODEL] --sampler SAMPLER
                                  [--start-time START_TIME]
                                  [--max-chain-encoding-k MAX_CHAIN_ENCODING_K]
                                  [--exact-prior] [--discard-ligand]
                                  [--discard-sdf-coords] [--detect-covalent]
                                  [--use-template] [--csv-path CSV_PATH]

    optional arguments:
      -h, --help            show this help message and exit
      --task TASK
      --sample-id SAMPLE_ID
      --template-id TEMPLATE_ID
      --cuda
      --model-checkpoint MODEL_CHECKPOINT
      --input-ligand INPUT_LIGAND
      --input-receptor INPUT_RECEPTOR
      --input-template INPUT_TEMPLATE
      --out-path OUT_PATH
      --n-samples N_SAMPLES
      --chunk-size CHUNK_SIZE
      --num-steps NUM_STEPS
      --latent-model LATENT_MODEL
      --sampler SAMPLER
      --start-time START_TIME
      --max-chain-encoding-k MAX_CHAIN_ENCODING_K
      --exact-prior
      --discard-ligand
      --discard-sdf-coords
      --detect-covalent
      --use-template
      --csv-path CSV_PATH

    """
    Path(output_path).mkdir(parents=True, exist_ok=True)
    jobs = []
    for job in job_data:
        command = ["apptainer", "run", "--nv", image, 
                   "neuralplexer-inference",
                   "--task=batched_structure_sampling"]
        command += ["--model-checkpoint", checkpoint]

        output = f"{output_path}/{job['name']}/"
        Path(output).mkdir(parents=True, exist_ok=True)
        command += ["--out-path", output]

        receptor_string = receptor(job['sequences']) 
        command += ["--input-receptor", f"\"{receptor_string}\""]

        ligand_string = ligand(job['ligands']) 
        command += ["--input-ligand", f"\"{ligand_string}\""]

        command += ["--sampler", job['parameters']['sampler']]
        command += ["--n-samples", job['parameters']['n-samples']]
        if cuda:
            command += ["--cuda"]
        jobs.append(command)
    return jobs


def find_next_item(lst, match_item):
    try:
        index = lst.index(match_item)
        if index + 1 < len(lst):
            return lst[index + 1]
        else:
            raise ValueError(f"Error: '{match_item}' is the last item in the" \
                             "list, no next item to return.")
    except ValueError as e:
        raise ValueError(f"Error: {str(e)}")


def run_jobs(jobs):
    # TODO paralelize this
    for job in jobs:
        cmd = " ".join(str(v) for v in job)
        print(cmd)
        subprocess.run(cmd, shell=True)


def main():
    """The main function for the neuralPlexer wrapper.

    Must be run within the neuralPlexer pixi evironment.
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('image', 
                        type=str, 
                        help="Path to neuralplexer image")
    parser.add_argument('model_checkpoint', 
                        type=str, 
                        help="Path to neuralplexer checkpoint")
    parser.add_argument('json_file', 
                        type=str, 
                        help="Path to the JSON file to be parsed")
    parser.add_argument('output_path',
                        type=str,
                        help="Where NP will write the results.")
    args = parser.parse_args()
    check_environment()

    job = parse_json(args.json_file)
    jobs = make_jobs(job, args.image, args.model_checkpoint, args.output_path)
    run_jobs(jobs)


if __name__ == "__main__":

    main()
