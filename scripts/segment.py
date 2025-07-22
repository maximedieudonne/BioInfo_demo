import imageio
import numpy as np
from skimage.filters import threshold_otsu
from skimage.io import imsave
import os
import sys

input_path = snakemake.input[0]
output_path = snakemake.output[0]

# Assure que le répertoire de sortie existe
output_dir = os.path.dirname(output_path)
os.makedirs(output_dir, exist_ok=True)

# Lecture de l'image
img = imageio.imread(input_path)
if img.ndim == 3:
    img = img.mean(axis=2)

# Segmentation simple par seuillage
thresh = threshold_otsu(img)
binary = img > thresh

# Sauvegarde de l'image segmentée
imsave(output_path, (binary * 255).astype(np.uint8))

# Vérification manuelle (debug)
if not os.path.exists(output_path):
    print(f"ERREUR : le fichier {output_path} n'a pas été créé.")
    sys.exit(1)
