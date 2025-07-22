configfile: "config.yaml"

rule all:
    input:
        expand("results/masks/{img}.png", img=[i.replace('.tif', '') for i in config["images"]]),
        expand("results/colored/{img}.png", img=[i.replace('.tif', '') for i in config["images"]])

rule segment:
    input: "data/raw/{img}.tif"
    output: "results/masks/{img}.png"
    conda: "envs/bioinfo_env.yaml"
    script: "scripts/segment.py"

rule colorize:
    input:
        img="data/raw/{img}.tif",
        mask="results/masks/{img}.png"
    output: "results/colored/{img}.png"
    conda: "envs/bioinfo_env.yaml"
    script: "scripts/colorize.py"

