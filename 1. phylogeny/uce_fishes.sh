# Find the probes in each alignement
phyluce_probe_run_multiple_lastzs_sqlite \
    --db results/fish/alignment_fish.sqlite \
    --output results/fish/alignment-fish-lastz \
    --scaffoldlist Ailurus_fulgens Amphiprion_ocellaris Aptenodytes_forsteri Ara_glaucogularis Arctocephalus_gazella Betta_splendens Bubo_scandiacus Callithrix_jacchus Canis_lupus_familiaris Capra_hircus Cavia_aperea Ceratotherium_simum_simum Cercocebus_lunulatus Cervus_elaphus_yarkandensis Cervus_nippon Chauna_torquata Chrysemys_picta Coleonyx_brevis Coturnix_japonica Cyanistes_caeruleus Cynoglossus_semilaevis Cyprinus_carpio Elephas_maximus Eublepharis_macularius Felis_catus Fukomys_damarensis Gallus_gallus Giraffa_camelopardalis Gyps_fulvus Hippopotamus_amphibius Homo_sapiens Hylobates_lar Larimichthys_crocea Larus_argentatus Larus_marinus Macaca_mulatta Mandrillus_leucophaeus Monodelphis_domestica Moschus_berezovskii Muntiacus_reevesi Mus_musculus Neovison_vison Odobenus_rosmarus Orcinus_orca Oryzias_latipes Panthera_pardus Panthera_tigris Pan_troglodytes Paralichthys_olivaceus Pelecanus_crispus Phoenicopterus_roseus Pithecia_pithecia Platalea_ajaja Pogona_vitticeps Procavia_capensis Pygoscelis_adeliae Rangifer_tarandus Rhea_pennata Rousettus_aegyptiacus Saimiri_boliviensis_boliviensis Salmo_salar Sarcophilus_harrisii Saxicola_maurus Sphaerodactylus_macrolepis Sus_scrofa Syngnathus_scovelli Taeniopygia_guttata Tapirus_indicus Thamnophis_sirtalis Tupaia_belangeri Turdus_merula Tursiops_truncatus Vicugna_pacos Vulpes_vulpes \
    --genome-base-path ./ \
    --probefile /home/lucie/PhyloMut/UCE/UCE_probes/Acanthomorphs-UCE-1Kv1-DUPE-SCREENED.fasta \
    --cores 12

## Extract the region in the ref with a flanking sequence of 1000 bp each side
phyluce_probe_slice_sequence_from_genomes \
    --lastz results/fish/alignment-fish-lastz \
    --conf genomes.conf \
    --flank 1000 \
    --name-pattern "Acanthomorphs-UCE-1Kv1-DUPE-SCREENED.fasta_v_{}.lastz.clean" \
    --output results/fish/alignment-fish-fasta

# Match back the contigs (fasta) to probes
phyluce_assembly_match_contigs_to_probes \
    --contigs results/fish/alignment-fish-fasta \
    --probes /home/lucie/PhyloMut/UCE/UCE_probes/Acanthomorphs-UCE-1Kv1-DUPE-SCREENED.fasta \
    --output results/fish/uce-search-results
mkdir -p results/fish/taxon-sets/all

## Create the data matrix to configurate the files
phyluce_assembly_get_match_counts \
    --locus-db results/fish/uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output results/fish/taxon-sets/all/all-taxa-incomplete.conf
mkdir results/fish/taxon-sets/all/log

# Get all the probes for all taxa in a fasta format all-taxa-incomplete.fasta
phyluce_assembly_get_fastas_from_match_counts \
    --contigs results/fish/alignment-fish-fasta \
    --locus-db results/fish/uce-search-results/probe.matches.sqlite \
    --match-count-output results/fish/taxon-sets/all/all-taxa-incomplete.conf \
    --output results/fish/taxon-sets/all/all-taxa-incomplete.fasta \
    --incomplete-matrix results/fish/taxon-sets/all/all-taxa-incomplete.incomplete \
    --log-path results/fish/taxon-sets/all/log

# With no trimming and mafft aligner, output one fasta for each probes
phyluce_align_seqcap_align \
    --fasta results/fish/taxon-sets/all/all-taxa-incomplete.fasta \
    --output results/fish/taxon-sets/all/mafft-nexus-internal-trimmed \
    --taxa 74 \
    --aligner mafft \
    --cores 6 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path results/fish/taxon-sets/all/log
    
# Get sequences by taxon
phyluce_assembly_explode_get_fastas_file \
    --input results/fish/taxon-sets/all/all-taxa-incomplete.fasta \
    --output results/fish/exploded-fastas \
    --by-taxon

# Clean the name of sequences
phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments results/fish/taxon-sets/all/mafft-nexus-internal-trimmed \
    --output results/fish/taxon-sets/all/mafft-nexus-internal-trimmed-clean \
    --input-format fasta
    --cores 6 \
    --log-path results/fish/taxon-sets/all/log/

# Create a 75% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments results/fish/taxon-sets/all/mafft-nexus-internal-trimmed-clean \
    --taxa 74 \
    --percent 0.75 \
    --output results/fish/taxon-sets/all/mafft-nexus-internal-trimmed-clean-75p \
    --cores 6 \
    --log-path results/fish/taxon-sets/all/log/

# Create a 80% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments results/fish/taxon-sets/all/mafft-nexus-internal-trimmed-clean \
    --taxa 74 \
    --percent 0.80 \
    --output results/fish/taxon-sets/all/mafft-nexus-internal-trimmed-clean-80p \
    --cores 6 \
    --log-path results/fish/taxon-sets/all/log/

# Create a 90% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments results/fish/taxon-sets/all/mafft-nexus-internal-trimmed-clean \
    --taxa 74 \
    --percent 0.90 \
    --output results/fish/taxon-sets/all/mafft-nexus-internal-trimmed-clean-90p \
    --cores 6 \
    --log-path results/fish/taxon-sets/all/log/


