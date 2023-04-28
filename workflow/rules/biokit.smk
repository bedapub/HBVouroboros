rule unmapped_sample_annotation:
    input:
        biokit_outdir = biokit_outdir
    output:
        output_file = "results/annot/phenoData.meta"
    run:
        biokit.biokit_unmapped_sample_annotation(biokit_outdir, output_file)
