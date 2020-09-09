from HBVouroboros import biokit

biokit_dir = config["biokit_dir"]
output_file = config["output_file"]

rule unmapped_sample_annotation:
    input:
        biokit_dir=biokit_dir
    output:
        output_file=output_file
    run:
        biokit.biokit_unmapped_sample_annotation(biokit_dir, output_file)
