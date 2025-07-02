# metadata read-in
meta_path = system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
meta = read.csv(meta_path)

# ANI/Abundance data read in
# sylph
sylph_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
se <- strainspy::read_sylph(sylph_path, meta)

# metaphlan
mp_path <- system.file("extdata", "metaphlan_merged.tsv.gz", package = "strainspy")
mp <- strainspy::read_metaphlan(mp_path, meta)

# sourmash
sm_path <- system.file("extdata", "example_sourmash.csv.gz", package = "strainspy")
sm <- strainspy::read_sourmash(sm_path, meta)

# taxonomy
tax_path <- system.file("extdata", "example_taxonomy.tsv.gz", package = "strainspy")
tax = read.csv(tax_path)

mp_tax_path <- system.file("extdata", "metaphlan_taxonomy.tsv.gz", package = "strainspy")
mp_tax = read.csv(mp_tax_path)


# Test a standard model 