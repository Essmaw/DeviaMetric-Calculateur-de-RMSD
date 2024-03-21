# DeviaMetric-Calculateur-de-RMSD

Bienvenue dans DeviaMetric !

Le RMSD (Root Mean Square Deviation) est une mesure couramment utilisée en bio-informatique pour évaluer la similarité structurelle entre deux ensembles d'atomes. Plus précisément, il quantifie la moyenne des distances entre les atomes de molécules superposées.

---

À propos de DeviaMetric 🌟

DeviaMetric est un outil conçu pour calculer le RMSD entre deux fichiers PDB, tout en offrant une visualisation interactive de l'alignement entre les structures moléculaires correspondantes. Mais, il calcule aussi le RMSD entre une structure de réference et une trajectoire de simulation !

## Fonctionnalités

- Calcul du RMSD entre deux structures PDB ou entre une structure de référence et une trajectoire de simulation.
- Choix des atomes sur lesquels baser l'alignement (tels que les carbones alpha, tous les atomes ou seulement le backbone).
- Visualisation en 3D à l'aide de NGLView.
- Interface conviviale pour charger et analyser les fichiers PDB.
- Possibilité d'installer les bibliothèques nécessaires via le fichier requirement.txt.
- Exécution du Jupyter Notebook via Binder pour une expérience interactive.

---

## Prérequis

Assurez-vous d'avoir les bibliothèques suivantes installées :

- Biopython
- NGLview
- Tkinter
- MDAnalysis
- Matplotlib
- Ipywidgets

Vous pouvez installer toutes ces bibliothèques en exécutant la commande suivante :

    ```bash
        pip install -r requirements.txt
    ```

---

## Comment exécuter

Depuis le répertoire principal de DeviaMetric, exécutez le fichier principal `interface.py` avec Python pour lancer l'application Tkinter :
    
    ```bash
        python3 src/Fixed-Structure/interface.py
    ```


Pour exécuter le Jupyter Notebook via Binder, [![Cliquer ici !](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Essmaw/DeviaMetric-Calculateur-de-RMSD/HEAD) 

---

Nous espérons que vous apprécierez l'utilisation de DeviaMetric pour vos besoins d'analyse en bio-informatique ! 

N'hésitez pas à nous faire part de vos retours et suggestions pour améliorer notre outil 🚀