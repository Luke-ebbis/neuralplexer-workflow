# Notes

## Running neural plexer

```bash

pixi apptainer run --nv results/dependencies/neuralplexer/neuralplexer.sif neuralplexer-inference --task=batched_structure_sampling --input-receptor programs/p2rank_2.4.1/prank_test/rcsb/1b8uA.pdb --sampler=langevin_simulated_annealing --cuda --model-checkpoint complex_structure_prediction.ckpt --out-path . --n-samples 16 --chunk-size 4 --num-steps=40 --input-ligand CID_5317570.sdf
```

 Update checkpoint
 python -m pytorch_lightning.utilities.upgrade_checkpoint --file results/dependencies/neuralplexer/data/neuralplexermodels_downstream_datasets_predictions/models/complex_structure_prediction.ckpt
