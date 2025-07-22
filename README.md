# Image Processing FAIR Pipeline

Ce projet propose un pipeline de traitement d’images de cellules au format `.tif`, construit avec **Snakemake** et conçu selon les principes **FAIR** (Findable, Accessible, Interoperable, Reusable).

Il inclut les étapes suivantes :
1. **Segmentation** : détection des cellules (premier plan / arrière-plan)
2. **Colorisation** : chaque cellule est affichée dans une couleur différente

---

## Structure du projet

cell_tracking_fair_project/
├── Snakefile # Définition du pipeline Snakemake
├── config.yaml # Configuration des images à traiter
├── envs/ # Environnements Conda pour chaque étape
│ └── bioinfo_env.yaml
├── data/ # Données d'entrée
│ ├── raw/ # Images TIFF d'origine
│ └── metadata.csv # Métadonnées des images
├── results/ # Résultats générés par le pipeline
│ ├── masks/
│ └── colored/
├── scripts/ # Scripts Python pour le traitement
│ ├── segment.py
│ └── colorize.py
├── docs/ # Rapport ou documentation complémentaire
├── LICENSE # Licence du projet
├── CITATION.cff # Fichier de citation
└── README.md # Ce fichier


---

## Lancer le pipeline

### Prérequis

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [WSL + Ubuntu](https://learn.microsoft.com/fr-fr/windows/wsl/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) installé dans un environnement minimal

### Installation

Clonez ce dépôt :

```bash
git clone https://github.com/votre-utilisateur/bioinfo_demo.git
cd bioinfo_demo
```

### Configuration

Modifiez le fichier config.yaml pour définir les images à traiter :

```yaml
images:
  - t001.tif
  - t002.tif
```
Les images correspondantes doivent se trouver dans data/raw/.

### Execution

snakemake --use-conda --cores 4 --latency-wait 20

## Environement Conda

Les dépendances sont gérées automatiquement par Snakemake à partir du fichier : envs/bioinfo_env.yaml

Snakemake créera les environnements requis automatiquement (avec --use-conda).

## Licence

Ce projet est distribué sous la licence MIT. Voir LICENSE.

## Citation

Si vous utilisez ce pipeline, merci de citer le projet :

```yaml
cff-version: 1.2.0
title: "Cell Tracking FAIR Pipeline"
authors:
  - family-names: Maxidieu
    given-names: <Prénom>
license: MIT
version: 1.0
```

Voir CITATION.cff.

## Contact

Pour toute question, contribution ou bug, contactez : maximedieudonne@protonmail.com