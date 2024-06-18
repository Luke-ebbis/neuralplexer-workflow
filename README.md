# NeuralPlexer workflow

Predicting the interaction between ligands and proteins with NeuralPlexer, a
deeplearning tool by [Qiao _et al._
2024](https://www.nature.com/articles/s42256-024-00792-z).

## Installation

The following commands install dependencies of the workflow in the _current_
directory within the `.pixi` folder. After installation, you cannot move the
folder without re-installling all the dependencies. 

```bash
curl -fsSL https://pixi.sh/install.sh | bash
# ... cd <this repo>
pixi install
```

## Usage

Each job has a `json` descriptor in the `data/` folder:

```json
[
  {
    "name": "2024-05-13_20:16",
    "parameters": {
      "sampler":"langevin_simulated_annealing",
      "n-samples": 10,
      "chunk-size": 1,
      "num-steps": 40
    },
    "ligands" : [
       {
        "molecule": {
          "sdf": "data/lig_ref.sdf",
          "count": 2
        }
      }
    ],
    "sequences": [
      {
        "proteinChain": {
          "sequence": "FGGGFGGGGGSGSGSGG",
          "count": 2
        }
      },
      {
        "proteinChain": {
          "sequence": "FGGSGSGSGG",
          "count": 1
        }
      }
    ]
  }
]
```

The ligands are included like this.

```
data
├── lig_ref.sdf
└── test.json
```

This will make a folder in `results/data/<jobname>` with each job in it. Right
now, the token limit seems to be 600 amino acid residues. At this point, the
graphical card runs out of memory.


* `pixi run make` runs the full workflow. You can supply arguments to
`snakemake` as needed, such as `--cores 10`, if your process needs 10 cores.

* `pixi run test` runs a dry run (`-np`) of the workflow. 

* `pixi run update_dag` updates the directed acyclic graph of the snakemake
workflow and moves it to the `resources` folder.

* `pixi run help` show the usage of the Snakefile


Usage in a HPC context is done with the

```
pixi run slurm
```

command, this will launch the hardware intensive snakemake tasks as SLURM jobs.

## About

When `pixi run make` is executed, the snakemake pipeline sets up neuralplexer
and predicts each listed complex in `data`. Additionally, summary statistics
are calculated and plotted as shown below. These summary statistics end up in
`results/analysis`. Statistics and analysis include:

* ramachandran plot.

* residue wise plddt.

* sasa surface Accessibility.

* Number of steric clashes.

* Identification of interacting residues using plip.

![](resources/pipeline.png)
