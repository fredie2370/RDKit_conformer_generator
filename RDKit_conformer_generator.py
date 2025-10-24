# analysis/utils/util.py
#MGR


import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem

def funcion(out, ligand, log_callback=None):

    def log(message):
        if log_callback:
            log_callback(message) 
        else:
            print(message)      

    out_dir = out
    ligando = ligand

    if not out_dir:
        log("Error: No se especificó el directorio de salida :(")
        return
    if not ligando:
        log("No se proporcionó un ligando :(")
        return

    log(f"Generando conformeros en: {out_dir}")
    if isinstance(ligand, tuple):
        name = ligando[0]
        smiles = ligando[1]
        log(f"Procesando ligando {name}, con SMILES {smiles}")
        run_rdkit(name, smiles, out_dir, log)
    else:
        log(f"Procesando archivo {ligando}")
        try:
            with open(ligando, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    smiles = parts[0]
                    name = "_".join(parts[1:])
                    run_rdkit(name, smiles, out_dir, log)
        except FileNotFoundError:
            log(f"Error. el archivo {ligando} no fue encontrado.")

def run_rdkit(name, smiles, out_dir, log):
    conformeros = 250
    max_final_conformers = 25
    rmsd_threshold = 1.0
    log(f"Generando pool de {conformeros} conformeros para {name}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log(f"No se pudo leer el SMILES para {name}")
        return 
    mol.SetProp("_Name", name)
    mol_h = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol_h,
                                    numConfs = conformeros,
                                    randomSeed=1,
                                    useExpTorsionAnglePrefs=True,
                                    useBasicKnowledge=True,
                                    numThreads=0)
    if len(cids) == 0:
        log(f"Advertencia: no se generaron conformeros para {name}. Saltando.")
        return 
    log(f"{len(cids)} conformeros generados. Optimizando...")
    try: 
        res = AllChem.MMFFOptimizeMoleculeConfs(mol_h, numThreads=0, maxIters=2000)
        log(f"{name} optimizado.")
    except Exception as e:
        log(f"     ADVERTENCIA: Falló la optimización MMFF94 para {name}: {e}")
    
    molecule_ouput_dir = os.path.join(out_dir, "Ligandos", name)
    os.makedirs(molecule_ouput_dir, exist_ok=True)
    
    output_filename = f"{name}_all_conformers.sdf"
    output_path = os.path.join(molecule_ouput_dir, output_filename)
    writer = Chem.SDWriter(output_path)
    for conf in mol_h.GetConformers():
        writer.write(mol_h, confId=conf.GetId())
    writer.close()
    log(f"     Éxito: {len(cids)} confórmeros guardados en {output_path}.")

    log(f"Generando {max_final_conformers} para {name}.")
    log(f"Threshold: {rmsd_threshold}.")

    energies = [tlp[1] for tlp in res if tlp [0] == 0]

    print(len(energies))

    #if len(energies) != mol_h.GetConformers():
    #    energies = [AllChem.MMFFGetMoleculeProperties(mol_h, "MMFF94s", confId=cid).GetMMFFEnergy() for cid in cids]
    sorted_confs = sorted(zip(cids, energies), key=lambda x: x[1])

    #print(sorted_confs)

    pruned_cids = []

    for cid, energy in sorted_confs:
        if len(pruned_cids) >= max_final_conformers:
            log(f"Limite de {max_final_conformers} confórmeros alcanzado.")
            break
        is_diverse =True

        for pruned_cid in pruned_cids:
            rmsd = AllChem.GetBestRMS(mol_h, mol_h, pruned_cid, cid)
            if rmsd < rmsd_threshold:
                is_diverse = False
                break
        if is_diverse:
            pruned_cids.append(cid)

    log(f"Filtrado completo {len(pruned_cids)} conformeros diversos seleccionados.")

    output_filename_diverse = f"{name}_diverse_conformers.sdf"
    output_path_diverse = os.path.join(molecule_ouput_dir, output_filename_diverse)

    writer_diverse = Chem.SDWriter(output_path_diverse)
    for cid in pruned_cids:
        writer_diverse.write(mol_h, confId=cid)
    writer_diverse.close()

    log(f"Exito: {len(pruned_cids)} conformeros guardados en {output_path_diverse}.")





if __name__ == "__main__":
    lig = ("nombre", "OC(=O)CCCCCOc1ccccc1CN(C(=O)c1ccc(cc1)c1ccco1)C1CC1")
    funcion("./", lig, log_callback=print)