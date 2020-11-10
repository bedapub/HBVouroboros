biokit_dir = config["biokit_dir"]

rule unmapped_sample_annotation:
    input:
        biokit_dir = biokit_dir
    output:
        output_file = "results/annot/phenoData.meta"
    run:
        biokit.biokit_unmapped_sample_annotation(biokit_dir, output_file)
