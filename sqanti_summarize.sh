#!/bin/bash

for assembler in flair stringtie.free stringtie.guided rb rattle
do
    echo analysis for ${assembler} >> sqanti_summary.txt

    wc -l *${assembler}.*classification.txt | awk '{print "transcripts assembled: " $1 -1}' >> sqanti_summary.txt # no. of transcripts (-1 for the header line)
    
    cut -f8 *${assembler}.*classification.txt | grep -v novel | sort | uniq | wc -l | \
    awk '{print "unique associated transcripts: " $1 -1}' >> sqanti_summary.txt #no. of unique associated transcripts
    
    grep 'full-splice_match' *${assembler}.*classification.txt | wc -l | awk '{print "FSM:\t" $1}' >> sqanti_summary.txt #for fsm
    grep 'full-splice_match' *${assembler}.*classification.txt | \
    cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM:\t" $1}' >> sqanti_summary.txt #for unique fsm

    grep 'full-splice_match' *${assembler}.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | \
    wc -l | awk '{print "RM:\t" $1}' >> sqanti_summary.txt #for rm
    grep 'full-splice_match' *${assembler}.*classification.txt | grep -e 'reference_match' -e 'mono-exon' | \
    cut -f8 | sort | uniq | wc -l | awk '{print "Unique RM:\t" $1}' >> sqanti_summary.txt #for unique rm

    grep 'incomplete-splice_match' *${assembler}.*classification.txt | wc -l | awk '{print "ISM:\t" $1}' >> sqanti_summary.txt #for ism
    grep 'incomplete-splice_match' *${assembler}.*classification.txt | \
    cut -f8 | sort | uniq | wc -l | awk '{print "Unique ISM:\t" $1}' >> sqanti_summary.txt #for unique ism

    grep -e 'full-splice_match' -e 'incomplete-splice_match' *${assembler}.*classification.txt | \
    cut -f8 | sort | uniq | wc -l | awk '{print "Unique FSM/ISM:\t" $1}' >> sqanti_summary.txt #for unique FSM/ISM

    grep 'novel_in_catalog' *${assembler}.*classification.txt | wc -l | awk '{print "NIC:\t" $1}' >> sqanti_summary.txt # for NIC
    grep 'novel_not_in_catalog' *${assembler}.*classification.txt | wc -l | awk '{print "NNC:\t" $1}' >> sqanti_summary.txt # for NNC
    grep 'antisense' *${assembler}.*classification.txt | wc -l | awk '{print "antisense:\t" $1}'>> sqanti_summary.txt # for antisense
    grep -w 'genic' *${assembler}.*classification.txt | wc -l | awk '{print "genic:\t" $1}' >> sqanti_summary.txt # for genic
    grep 'intergenic' *${assembler}.*classification.txt | wc -l | awk '{print "intergenic:\t" $1}' >> sqanti_summary.txt # for intergenic
    grep 'fusion' *${assembler}.*classification.txt | wc -l | awk '{print "fusion:\t" $1}' >> sqanti_summary.txt # for fusion

done
