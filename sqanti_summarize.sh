#!/bin/bash

# Copy this script to the assembly that you want to analyze in 03.sqanti
echo 'analysis for flair' > sqanti_summary.txt

wc -l *flair.*classification.txt | awk '{print "transcripts assembled: " $1 -1}' >> sqanti_summary.txt # no. of transcripts (-1 for the header line)

grep 'full-splice_match' *flair.*classification.txt | wc -l | awk '{print "FSM:\t" $1}' >> sqanti_summary.txt #for fsm
grep 'full-splice_match' *flair.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM:\t" $1}' >> sqanti_summary.txt #for unique fsm

grep 'full-splice_match' *flair.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | \
wc -l | awk '{print "RM:\t" $1}' >> sqanti_summary.txt #for rm
grep 'full-splice_match' *flair.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique RM:\t" $1}' >> sqanti_summary.txt #for unique rm

grep 'incomplete-splice_match' *flair.*classification.txt | wc -l | awk '{print "ISM:\t" $1}' >> sqanti_summary.txt #for ism
grep 'incomplete-splice_match' *flair.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique ISM:\t" $1}' >> sqanti_summary.txt #for unique ism

grep -e 'full-splice_match' -e 'incomplete-splice_match' *flair.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM/ISM:\t" $1}' >> sqanti_summary.txt #for unique FSM/ISM

grep 'novel_in_catalog' *flair.*classification.txt | wc -l | awk '{print "NIC:\t" $1}' >> sqanti_summary.txt # for NIC
grep 'novel_not_in_catalog' *flair.*classification.txt | wc -l | awk '{print "NNC:\t" $1}' >> sqanti_summary.txt # for NNC
grep 'antisense' *flair.*classification.txt | wc -l | awk '{print "antisense:\t" $1}'>> sqanti_summary.txt # for antisense
grep 'genic' *flair.*classification.txt | wc -l | awk '{print "genic:\t" $1}' >> sqanti_summary.txt # for genic
grep 'intergenic' *flair.*classification.txt | wc -l | awk '{print "intergenic:\t" $1}' >> sqanti_summary.txt # for intergenic
grep 'fusion' *flair.*classification.txt | wc -l | awk '{print "fusion:\t" $1}' >> sqanti_summary.txt # for fusion


echo 'analysis for stringtie.guided' >> sqanti_summary.txt

wc -l *stringtie.guided.*classification.txt | awk '{print "transcripts assembled: " $1 -1}' >> sqanti_summary.txt # no. of transcripts (-1 for the header line)

grep 'full-splice_match' *stringtie.guided.*classification.txt | wc -l | awk '{print "FSM:\t" $1}' >> sqanti_summary.txt #for fsm
grep 'full-splice_match' *stringtie.guided.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM:\t" $1}' >> sqanti_summary.txt #for unique fsm

grep 'full-splice_match' *stringtie.guided.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | \
wc -l | awk '{print "RM:\t" $1}' >> sqanti_summary.txt #for rm
grep 'full-splice_match' *stringtie.guided.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique RM:\t" $1}' >> sqanti_summary.txt #for unique rm

grep 'incomplete-splice_match' *stringtie.guided.*classification.txt | wc -l | awk '{print "ISM:\t" $1}' >> sqanti_summary.txt #for ism
grep 'incomplete-splice_match' *stringtie.guided.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique ISM:\t" $1}' >> sqanti_summary.txt #for unique ism

grep -e 'full-splice_match' -e 'incomplete-splice_match' *stringtie.guided.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM/ISM:\t" $1}' >> sqanti_summary.txt #for unique rm\ism

grep 'novel_in_catalog' *stringtie.guided.*classification.txt | wc -l | awk '{print "NIC:\t" $1}' >> sqanti_summary.txt # for NIC
grep 'novel_not_in_catalog' *stringtie.guided.*classification.txt | wc -l | awk '{print "NNC:\t" $1}' >> sqanti_summary.txt # for NNC
grep 'antisense' *stringtie.guided.*classification.txt | wc -l | awk '{print "antisense:\t" $1}'>> sqanti_summary.txt # for antisense
grep 'genic' *stringtie.guided.*classification.txt | wc -l | awk '{print "genic:\t" $1}' >> sqanti_summary.txt # for genic
grep 'intergenic' *stringtie.guided.*classification.txt | wc -l | awk '{print "intergenic:\t" $1}' >> sqanti_summary.txt # for intergenic
grep 'fusion' *stringtie.guided.*classification.txt | wc -l | awk '{print "fusion:\t" $1}' >> sqanti_summary.txt # for fusion


echo 'analysis for stringtie.free' >> sqanti_summary.txt

wc -l *stringtie.free.*classification.txt | awk '{print "transcripts assembled: " $1 -1}' >> sqanti_summary.txt # no. of transcripts (-1 for the header line)
grep 'full-splice_match' *stringtie.free.*classification.txt | wc -l | awk '{print "FSM:\t" $1}' >> sqanti_summary.txt #for fsm

grep 'full-splice_match' *stringtie.free.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM:\t" $1}' >> sqanti_summary.txt #for unique fsm
grep 'full-splice_match' *stringtie.free.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | \
wc -l | awk '{print "RM:\t" $1}' >> sqanti_summary.txt #for rm
grep 'full-splice_match' *stringtie.free.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique RM:\t" $1}' >> sqanti_summary.txt #for unique rm

grep 'incomplete-splice_match' *stringtie.free.*classification.txt | wc -l | awk '{print "ISM:\t" $1}' >> sqanti_summary.txt #for ism
grep 'incomplete-splice_match' *stringtie.free.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique ISM:\t" $1}' >> sqanti_summary.txt #for unique ism

grep -e 'full-splice_match' -e 'incomplete-splice_match' *stringtie.free.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM/ISM:\t" $1}' >> sqanti_summary.txt #for unique FSM/ISM

grep 'novel_in_catalog' *stringtie.free.*classification.txt | wc -l | awk '{print "NIC:\t" $1}' >> sqanti_summary.txt # for NIC
grep 'novel_not_in_catalog' *stringtie.free.*classification.txt | wc -l | awk '{print "NNC:\t" $1}' >> sqanti_summary.txt # for NNC
grep 'antisense' *stringtie.free.*classification.txt | wc -l | awk '{print "antisense:\t" $1}'>> sqanti_summary.txt # for antisense
grep 'genic' *stringtie.free.*classification.txt | wc -l | awk '{print "genic:\t" $1}' >> sqanti_summary.txt # for genic
grep 'intergenic' *stringtie.free.*classification.txt | wc -l | awk '{print "intergenic:\t" $1}' >> sqanti_summary.txt # for intergenic
grep 'fusion' *stringtie.free.*classification.txt | wc -l | awk '{print "fusion:\t" $1}' >> sqanti_summary.txt # for fusion


echo 'analysis for rb' >> sqanti_summary.txt

wc -l *rb.*classification.txt | awk '{print "transcripts assembled: " $1 -1}' >> sqanti_summary.txt # no. of transcripts (-1 for the header line)

grep 'full-splice_match' *rb.*classification.txt | wc -l | awk '{print "FSM:\t" $1}' >> sqanti_summary.txt #for fsm
grep 'full-splice_match' *rb.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM:\t" $1}' >> sqanti_summary.txt #for unique fsm

grep 'full-splice_match' *rb.*classification.txt | \
grep -e 'reference_match' -e 'mono-exon' | wc -l | awk '{print "RM:\t" $1}' >> sqanti_summary.txt #for rm
grep 'full-splice_match' *rb.*classification.txt | \
grep -e 'reference_match' -e 'mono-exon' | cut -f8 | sort | uniq | wc -l | awk '{print "Unique RM:\t" $1}' >> sqanti_summary.txt #for unique rm

grep 'incomplete-splice_match' *rb.*classification.txt | wc -l | awk '{print "ISM:\t" $1}' >> sqanti_summary.txt #for ism
grep 'incomplete-splice_match' *rb.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique ISM:\t" $1}' >> sqanti_summary.txt #for unique ism

grep -e 'full-splice_match' -e 'incomplete-splice_match' *rb.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM/ISM:\t" $1}' >> sqanti_summary.txt #for unique FSM/ISM

grep 'novel_in_catalog' *rb.*classification.txt | wc -l | awk '{print "NIC:\t" $1}' >> sqanti_summary.txt # for NIC
grep 'novel_not_in_catalog' *rb.*classification.txt | wc -l | awk '{print "NNC:\t" $1}' >> sqanti_summary.txt # for NNC
grep 'antisense' *rb.*classification.txt | wc -l | awk '{print "antisense:\t" $1}'>> sqanti_summary.txt # for antisense
grep 'genic' *rb.*classification.txt | wc -l | awk '{print "genic:\t" $1}' >> sqanti_summary.txt # for genic
grep 'intergenic' *rb.*classification.txt | wc -l | awk '{print "intergenic:\t" $1}' >> sqanti_summary.txt # for intergenic
grep 'fusion' *rb.*classification.txt | wc -l | awk '{print "fusion:\t" $1}' >> sqanti_summary.txt # for fusion


echo 'analysis for rattle' >> sqanti_summary.txt

wc -l *rattle.*classification.txt | awk '{print "transcripts assembled: " $1 -1}' >> sqanti_summary.txt # no. of transcripts (-1 for the header line)

grep 'full-splice_match' *rattle.*classification.txt | wc -l | awk '{print "FSM:\t" $1}' >> sqanti_summary.txt #for fsm
grep 'full-splice_match' *rattle.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM:\t" $1}' >> sqanti_summary.txt #for unique fsm
grep 'full-splice_match' *rattle.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | wc -l | awk '{print "RM:\t" $1}' >> sqanti_summary.txt #for rm
grep 'full-splice_match' *rattle.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique RM:\t" $1}' >> sqanti_summary.txt #for unique rm

grep 'incomplete-splice_match' *rattle.*classification.txt | wc -l | awk '{print "ISM:\t" $1}' >> sqanti_summary.txt #for ism
grep 'incomplete-splice_match' *rattle.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique ISM:\t" $1}' >> sqanti_summary.txt #for unique ism

grep -e 'full-splice_match' -e 'incomplete-splice_match' *rattle.*classification.txt | \
cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM/ISM:\t" $1}' >> sqanti_summary.txt #for unique FSM/ISM

grep 'novel_in_catalog' *rattle.*classification.txt | wc -l | awk '{print "NIC:\t" $1}' >> sqanti_summary.txt # for NIC
grep 'novel_not_in_catalog' *rattle.*classification.txt | wc -l | awk '{print "NNC:\t" $1}' >> sqanti_summary.txt # for NNC
grep 'antisense' *rattle.*classification.txt | wc -l | awk '{print "antisense:\t" $1}'>> sqanti_summary.txt # for antisense
grep 'genic' *rattle.*classification.txt | wc -l | awk '{print "genic:\t" $1}' >> sqanti_summary.txt # for genic
grep 'intergenic' *rattle.*classification.txt | wc -l | awk '{print "intergenic:\t" $1}' >> sqanti_summary.txt # for intergenic
grep 'fusion' *rattle.*classification.txt | wc -l | awk '{print "fusion:\t" $1}' >> sqanti_summary.txt # for fusion