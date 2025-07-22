import imageio
import numpy as np
from skimage.measure import label
from skimage.color import label2rgb
from skimage.io import imsave
import os
import sys

# Entrées et sorties
img_path = snakemake.input.img
mask_path = snakemake.input.mask
output_path = snakemake.output[0]

# Crée le dossier de sortie si besoin
output_dir = os.path.dirname(output_path)
os.makedirs(output_dir, exist_ok=True)

# Lecture du masque binaire
mask = imageio.imread(mask_path)
binary_mask = mask > 0

# Labellisation des régions
labels = label(binary_mask)

# Coloration des régions avec des couleurs aléatoires
colored = label2rgb(labels, bg_label=0, kind='overlay')

# Sauvegarde de l’image coloriée
imsave(output_path, (colored * 255).astype(np.uint8))

# Vérification que le fichier a bien été écrit
if not os.path.exists(output_path):
    print(f"ERREUR : le fichier {output_path} n'a pas été créé.")
    sys.exit(1)
