#!/usr/bin/python
# This script reads the TF database files and outputs the unique genes.

lines = open("inst/extdata/c3.tft.v4.0.entrez.gmt","r").read().split("\n")
genes = []
for line in lines:
	fields = line.split("\t")
	genes.extend(fields[2:])
genes = set(genes)
for gene in genes:
	print gene
