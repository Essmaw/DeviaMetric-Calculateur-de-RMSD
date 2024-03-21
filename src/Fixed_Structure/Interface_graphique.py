"""RMSD Calculator Tool.

This Tkinter application calculates the Root-Mean-Square Deviation (RMSD) between two trajectories or PDB files.
It allows users to select the files and displays the average RMSD in Angstroms. 


Usage:
======
From RMSD_Calculator/ repository, run:
    python run src/main.py

"""

__author__ = "Essmay TOUAMI"
__contact__ = "essmay.touami@etu.u-paris.fr"
__copyright__ = "Universite Paris-Cite"
__date__ = "2024-03-18"

##################
# import modules #
##################
import os
import warnings
import urllib.request
import tkinter as tk
from tkinter import ttk
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, Superimposer
import pyperclip
# pour ignorer les warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.simplefilter('ignore', BiopythonWarning)

################
# define class #
################
class RMSDCalculatorApp:
    def __init__(self, master):
        self.master = master
        master.title("DeviaMetric-Structure fixe")

        # Créer les widgets pour les variables
        self.pdb_id_1_label = ttk.Label(master, text="ID PDB 1:")
        self.pdb_id_1 = ttk.Entry(master)

        self.pdb_id_2_label = ttk.Label(master, text="ID PDB 2:")
        self.pdb_id_2 = ttk.Entry(master)

        self.option_ali_label = ttk.Label(master, text="Type d'atomes pour l'alignement:")
        self.option_ali = ttk.Combobox(master, values=['CA', 'all', 'backbone'])

        # Créer le bouton pour calculer et afficher le RMSD
        self.calculate_button = ttk.Button(master, text="Calculer et Afficher RMSD", command=self.calculate_and_display_rmsd)

        # Créer la zone de texte pour afficher le RMSD
        self.output_text = tk.Text(master, height=5, width=60)

        # Créer le bouton pour copier le RMSD dans le presse-papiers
        self.copy_button = ttk.Button(master, text="Copier RMSD", command=self.copy_rmsd_to_clipboard)

        # Placer les widgets dans la grille
        self.pdb_id_1_label.grid(row=0, column=0, padx=5, pady=5)
        self.pdb_id_1.grid(row=0, column=1, padx=5, pady=5)

        self.pdb_id_2_label.grid(row=1, column=0, padx=5, pady=5)
        self.pdb_id_2.grid(row=1, column=1, padx=5, pady=5)

        self.option_ali_label.grid(row=2, column=0, padx=5, pady=5)
        self.option_ali.grid(row=2, column=1, padx=5, pady=5)

        self.calculate_button.grid(row=3, column=0, columnspan=2, padx=5, pady=5)
        self.output_text.grid(row=4, column=0, columnspan=2, padx=5, pady=5)
        self.copy_button.grid(row=5, column=0, columnspan=2, padx=5, pady=5)

    def get_structure_from_pdb_id(self,pdb_id):
        pdb_id = pdb_id.get()
        print(f"Téléchargement du fichier PDB de l'ID {pdb_id} en cours...")
        # Construction de l'URL PDB
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'

        try:
            # Télécharger le fichier PDB à partir de l'URL
            with urllib.request.urlopen(url) as response:
                pdb_content = response.read()

            # Définir le chemin local du fichier PDB
            pdb_file_path = f'{pdb_id.lower()}.pdb'

            # Écrire le contenu du fichier PDB dans le fichier local
            with open(pdb_file_path, 'wb') as pdb_file:
                pdb_file.write(pdb_content)
            print(f"Le fichier PDB de l'ID {pdb_id} a été téléchargé avec succès.")
            print("")
            return os.path.abspath(pdb_file_path)
        
        except urllib.error.HTTPError as e:
            self.output_text.insert(tk.END,f"Erreur lors de la récupération de la structure pour l'ID PDB {pdb_id}: {e}")
            return None


    def align_pdb_files(self, file1, file2, align_type):
        # Initialiser le parseur PDB
        parser = PDBParser()

        # Analyser les fichiers PDB et obtenir les structures
        self.structure1 = parser.get_structure('structure1', file1)
        self.structure2 = parser.get_structure('structure2', file2)

        # Accumuler les atomes pour l'alignement
        atoms1 = []
        atoms2 = []

        # Parcourir les résidus des deux structures
        for model1, model2 in zip(self.structure1, self.structure2):
            for chain1, chain2 in zip(model1, model2):
                for residue1, residue2 in zip(chain1, chain2):
                    # Sélectionner les atomes en fonction du type d'alignement
                    if align_type == "CA":
                        # Sélectionner les carbones alpha s'ils existent
                        if 'CA' in residue1 and 'CA' in residue2:
                            atoms1.append(residue1['CA'])
                            atoms2.append(residue2['CA'])
                        
                    elif align_type == "all":
                        # Sélectionner tous les atomes des résidus
                        atoms1.extend(residue1.get_atoms())
                        atoms2.extend(residue2.get_atoms())
                    elif align_type == "backbone":
                        # Sélectionner les atomes du backbone (N, CA, C, O)
                        backbone_atoms1 = [atom for atom in residue1 if atom.id in ['N', 'CA', 'C', 'O']]
                        backbone_atoms2 = [atom for atom in residue2 if atom.id in ['N', 'CA', 'C', 'O']]
                        # Ajouter les atomes du backbone à la liste
                        atoms1.extend(backbone_atoms1)
                        atoms2.extend(backbone_atoms2)

        # Vérifier que les listes d'atomes ont la même longueur
        min_length = min(len(atoms1), len(atoms2))
        atoms1 = atoms1[:min_length]
        atoms2 = atoms2[:min_length]

        # Créer l'objet Superimposer et définir les atomes
        sup = Superimposer()
        sup.set_atoms(atoms1, atoms2)

        # Appliquer la transformation aux atomes de la structure 2
        sup.apply(self.structure2.get_atoms())

        return self.structure1, self.structure2

    def calculate_rmsd(self, structure1, structure2, align_type):
        atoms1 = []
        atoms2 = []

        for model1, model2 in zip(structure1, structure2):
            for chain1, chain2 in zip(model1, model2):
                for residue1, residue2 in zip(chain1, chain2):
                    if align_type == "CA":
                        try:
                            atoms1.append(residue1['CA'])
                            atoms2.append(residue2['CA'])
                        except KeyError:
                            pass
                    elif align_type == "all":
                        atoms1.extend(residue1.get_atoms())
                        atoms2.extend(residue2.get_atoms())
                    elif align_type == "backbone":
                        # Sélectionner les atomes du backbone (N, CA, C, O)
                        backbone_atoms1 = [atom for atom in residue1 if atom.id in ['N', 'CA', 'C', 'O']]
                        backbone_atoms2 = [atom for atom in residue2 if atom.id in ['N', 'CA', 'C', 'O']]
                        atoms1.extend(backbone_atoms1)
                        atoms2.extend(backbone_atoms2)

        total_atoms = len(atoms1)
        rmsd = sum([(atom1 - atom2) ** 2 for atom1, atom2 in zip(atoms1, atoms2)]) / total_atoms
        rmsd = rmsd ** 0.5

        return rmsd
    
    def calculate_and_display_rmsd(self):
        # Ici, vous calculeriez le RMSD entre les structures et stockeriez la valeur dans la variable rmsd
        # on récupère les fichiers pdb téléchargés
        fichier1_pdb =  self.get_structure_from_pdb_id(self.pdb_id_1)
        fichier2_pdb =self.get_structure_from_pdb_id(self.pdb_id_2)

        # on aligne les structures selon le type d'atome selectionné
        structure1, structure2 = self.align_pdb_files(fichier1_pdb, fichier2_pdb, self.option_ali.get())
        
        # on calcule le rmsd
        rmsd = self.calculate_rmsd(structure1, structure2, self.option_ali.get())

        # Effacer le contenu précédent de la zone de texte
        self.output_text.delete(1.0, tk.END)
        
        # Affichage formaté du RMSD dans la zone de texte avec HTML
        output = f"Le RMSD entre la structure {self.pdb_id_1.get()} et {self.pdb_id_2.get()} est égale " \
                 f"à {rmsd:.4f} Å, en prenant en compte les atomes de type {self.option_ali.get()}."
        self.output_text.insert(tk.END, output)


    def copy_rmsd_to_clipboard(self):
        # Récupérer la valeur du RMSD depuis la zone de texte
        rmsd_value = self.output_text.get("1.0", tk.END)
        
        # Extraire le RMSD de la chaîne de texte
        rmsd = rmsd_value.split()[7]
        
        # Copier le RMSD dans le presse-papiers
        pyperclip.copy(rmsd)
        print("La valeur du RMSD a été copiée dans le presse-papiers !")


################
# main program #
################
if __name__ == "__main__":
    root = tk.Tk()
    root.resizable(False, False)
    app = RMSDCalculatorApp(root)
    root.mainloop()
