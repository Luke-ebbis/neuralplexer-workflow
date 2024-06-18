from Bio.PDB import PDBParser, PDBIO
from rdkit import Chem
import pathlib
from rdkit.Chem.rdmolfiles import SDMolSupplier

def merge_records(folder: str,
                  output: str):
    """Merge sdf records with PDB records.

    Parameters
    ----------
    folder: str
        The folder that holds sdf files and pdb files. Ligand files must be
        named `lig_{number}.sdf`. There must be a single file named
        `prot_all.pdb`.
    output: str
        Where the merged files need to be written.

    Code snippet provided by Hanna.
    """
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    parser = PDBParser()
    all_prot_struct = parser.get_structure("All_proteins",
                                           f"{folder}/prot_all.pdb")
    models = all_prot_struct.get_models()
    idx = 0
    for model in models:
        mol_supplier = SDMolSupplier(f"{folder}/lig_{idx}.sdf")
        for mol in mol_supplier:
            Chem.MolToPDBFile(mol, f"{folder}/lig_{idx}.pdb")
        lig_struc = parser.get_structure(f"Ligand{idx}",
                                         f"{folder}/lig_{idx}.pdb")
        chain = [chain for chain in lig_struc.get_chains()][0]
        chain.detach_parent()
        model.add(chain)
        
        io=PDBIO()
        io.set_structure(model)
        io.save(f"{output}/prot_lig_{idx}.pdb")
        idx += 1

def main():
    folder = snakemake.input[0]
    output = snakemake.output[0]
    merge_records(folder, output)
    

if __name__ == "__main__":

    main()
