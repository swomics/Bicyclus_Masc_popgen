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
