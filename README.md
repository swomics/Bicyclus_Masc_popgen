# Bicyclus_ambig_popgen
A record of finalised commands and scripts used

1. Haplotypes_unaligned_nt.fasta is pulled from Arjen's spreadsheet data. The sequences consist of the raw exon 8 and exon 9 sequences seperated by an arbitrary spacer consisting of 6 N's. African lab reared samples (F1's) are collapsed into singletons since they artificially increase the number of redundant alleles. Similarly singleton haplotypes from the resulting Leiden and Liverpool stocks are included as 'pseudoindividuals' to depict their position within the overall haplotype network.

2. Haplotypes_unaligned_trimmed_nt.fasta removes two 5' nucleotides from exon 8 and two 5' nucleotides from exon 9 to enable treatment as coding sequence.

3. Align with MACSE 

`java -jar ~/bin/macse_v2.03.jar -prog alignSequences -seq Haplotypes_unaligned_trimmed_nt.fasta`

4. Remove the HVR region of the alignment for recoding (AA position 23 to 54 or the fully conserved I residue to the fully conserved E residue)

5. Recode nucleotide and amino acid HVR usinf PICS-ord.

nucleotide params: `"-m zeta --free-end-gaps -o dist-c:-"`

amino acid params: `"-m aazeta --free-end-gaps -o dist-c:-"`

`R-3.3.3/bin/Rscript ./picsord_nt.R Haplotypes_MACSE_HVR_NT.fas > Haplotypes_MACSE_HVR_NT.phy`

`R-3.3.3/bin/Rscript ./picsord_aa.R Haplotypes_MACSE_HVR_AA.fas > Haplotypes_MACSE_HVR_AA.phy`

6. Convert and concatenate the 5prime and 3prime flanking alignments

nt

`AMAS.py convert -i 5prime_Haplotypes_unaligned_trimmed_nt_NT.fas -f fasta -d dna -u phylip`
 
`AMAS.py convert -i 3prime_Haplotypes_unaligned_trimmed_nt_NT.fas -f fasta -d dna -u phylip`
 
`python3 Phy_cat.py 5prime_Haplotypes_unaligned_trimmed_nt_NT.fas-out.phy Haplotypes_MACSE_HVR_NT.phy 3prime_Haplotypes_unaligned_trimmed_nt_NT.fas-out.phy`

aa

`AMAS.py convert -i 5prime_Haplotypes_unaligned_trimmed_aa_AA.fas -f fasta -d dna -u phylip`
 
`AMAS.py convert -i 3prime_Haplotypes_unaligned_trimmed_aa_AA.fas -f fasta -d dna -u phylip`
 
 concatenating this way caused problems as the partitions are too small
`### python3 Phy_cat.py 5prime_Haplotypes_unaligned_trimmed_aa_AA.fas-out.phy Haplotypes_MACSE_HVR_AA.phy 3prime_Haplotypes_unaligned_trimmed_aa_AA.fas-out.phy`

`python3 Phy_cat.py Haplotypes_MACSE_HVR_AA.phy 5prime_Haplotypes_unaligned_trimmed_aa_AA.fas-out.phy 3prime_Haplotypes_unaligned_trimmed_aa_AA.fas-out.phy`

`python3 Phy_cat.py Haplotypes_MACSE_HVR_NT.phy 5prime_Haplotypes_unaligned_trimmed_nt_NT.fas-out.phy 3prime_Haplotypes_unaligned_trimmed_nt_NT.fas-out.phy`

7aa. RaxML distances, raxml complains that only 19 states are present (no W), it is unclear how to resolve this as we cannot enlarge the partition.

`raxmlHPC-PTHREADS-SSE3 -q aa_part.txt -s superalgn_aa_rearr.fasta -K ORDERED -m PROTGAMMAAUTO -n T3 -p 12345 -f x`

`raxmlHPC-PTHREADS-SSE3 -q aa_part.txt -s superalgn_aa_rearr.fasta -K GTR -m PROTGAMMAAUTO -n T4 -p 12345 -f x`

`raxmlHPC-PTHREADS-SSE3 -q aa_part.txt -s superalgn_aa_rearr.fasta -K MK -m PROTGAMMAAUTO -n T4 -p 12345 -f x`


8aa. Examine the correlations between amino acid distance with no HVR vs raxml combined distances with the 3 different multi-state models and keep the model with the best correlation. We expect that the HVR contains phylogenetic signal that may differ from the flanks, however we make an assumption that the model that gives the best correlation is the most appropriate one for the HVR. We find GTR and MK models nearly identical, but GTR is the best (cor=0.517021, p-value < 2.2e-16), so we proceed with distances calculated with this model. Additionally, the GTR multi-state model also gives the highest final log-liklihood value.

7nt. RaxML distances

`raxmlHPC-PTHREADS-SSE3 -q nt_part.txt -s superalgn_nt_rearr.fasta -K GTR -m GTRGAMMA -n DNAT1 -p 12345 -f x`
