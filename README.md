# RDKit_conformer_generator

Generates an SDF conformer ensemble from SMILES. This script creates an initial pool of 250 conformers (using RDKit ETKDG), optimizes them (MMFF94), and then applies a greedy diversity filter. The final ensemble contains a subset of conformers selected for both their low energy and geometric diversity (RMSD threshold).

This aproach is based on: Andrew T. McNutt, Fatimah Bisiriyu, Sophia Song, Ananya Vyas, Geoffrey R. Hutchison, and David Ryan Koes. Journal of Chemical Information and Modeling 2023 63 (21), 6598-6607 DOI: 10.1021/acs.jcim.3c01245

## Usage

 In case you have de smiles in a .smi file:
    funcion("./", "/home/nevermore/Documentos/aplicacion_de_manu/Proyecto_bioinformático/scripts/pruebas/pru.smi" )

 In case you have a single smmiles 
    lig = ("nombre", "OC(=O)CCCCCOc1ccccc1CN(C(=O)c1ccc(cc1)c1ccco1)C1CC1")
    funcion("./", lig, log_callback=print)

The script will generate two files:
        * `*_conformers_all.sdf`: The full pool of 250 conformers.
        * `*_conformers_diverse.sdf`: The filtered set of diverse, low-energy conformers, ready for docking.


## License

This project is licensed under the GNU License. See the `LICENSE` file for details
