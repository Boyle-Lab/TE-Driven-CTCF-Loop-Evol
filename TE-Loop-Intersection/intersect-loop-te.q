-- RAD21 ChIA-pet data

CREATE TABLE ctcf_repeats_cp ROW FORMAT DELIMITED FIELDS TERMINATED BY '\t' STORED AS RCFile AS SELECT cs.*, cp.id AS `id_cp`, cp.start1, cp.end1, cp.start2, cp.end2, cp.iab, cp.fdr FROM ctcf_repeats cs INNER JOIN chia_pet cp ON (cp.cell = cs.cell) WHERE cp.chrom1 = cs.chrom AND cp.chrom2 = cs.chrom AND cp.factor='RAD21' AND ((cs.peak_ctcf BETWEEN cp.start1 AND cp.end1) OR (cs.peak_ctcf BETWEEN cp.start2 AND cp.end2));

CREATE TABLE ctcf_repeats_cp_mt ROW FORMAT DELIMITED FIELDS TERMINATED BY '\t' STORED AS RCFile AS SELECT cs.*, mt.id AS `id_mt`, mt.chromStart AS `chromStart_mt`, mt.chromEnd AS `chromEnd_mt`, mt.strand AS `strand_mt`, mt.score AS `score_mt` FROM ctcf_repeats_cp cs INNER JOIN motifs mt ON (mt.species = cs.species) WHERE cs.factor='CTCF' AND mt.factor = cs.factor AND mt.chrom = cs.chrom AND ( (mt.chromStart + mt.chromEnd)/2 BETWEEN cs.peak_ctcf-50 AND cs.peak_ctcf+50);

CREATE TABLE ctcf_cp ROW FORMAT DELIMITED FIELDS TERMINATED BY '\t' STORED AS RCFile AS SELECT cp.id AS `id_cp`, cp.chrom1 AS `chrom`, cp.start1, cp.end1, cp.start2, cp.end2, cp.iab, cp.fdr, cp.cell AS `cell`, cp.species AS `species`, cs.id AS `id_ctcf`, cs.chromstart AS `chromStart_ctcf`, cs.chromEnd AS `chromEnd_ctcf`, cs.signalValue AS `signalValue_ctcf`, cs.peak AS `peak_ctcf` FROM chia_pet cp INNER JOIN chip_seq cs ON (cs.cell = cp.cell) WHERE  cp.chrom1 = cs.chrom AND cp.chrom2 = cs.chrom AND cp.factor = 'RAD21' AND cs.factor='CTCF' AND ((cs.peak BETWEEN cp.start1 AND cp.end1) OR (cs.peak BETWEEN cp.start2 AND cp.end2));

CREATE TABLE ctcf_cp_mt ROW FORMAT DELIMITED FIELDS TERMINATED BY '\t' STORED AS RCFile AS SELECT cp.*, mt.id AS `id_mt`, mt.chromstart AS `start_mt`, mt.chromend AS `end_mt`, mt.score AS `score_mt`, mt.strand AS `strand_mt` FROM ctcf_cp cp INNER JOIN motifs mt ON (mt.species = cp.species) WHERE mt.factor = 'CTCF' AND mt.chrom = cp.chrom AND ( (mt.chromStart + mt.chromEnd)/2 BETWEEN cp.peak_ctcf-50 AND cp.peak_ctcf+50);

INSERT OVERWRITE DIRECTORY '/user/adadiehl/ctcf_te_mt_intersection' ROW FORMAT DELIMITED FIELDS TERMINATED BY'\t' STORED AS TEXTFILE SELECT * FROM ctcf_repeats_cp_mt;
INSERT OVERWRITE DIRECTORY '/user/adadiehl/ctcf_mt_intersection' ROW FORMAT DELIMITED FIELDS TERMINATED BY'\t' STORED AS TEXTFILE SELECT * FROM ctcf_cp_mt;
INSERT OVERWRITE DIRECTORY '/user/adadiehl/ctcf_te_intersection' ROW FORMAT DELIMITED FIELDS TERMINATED BY'\t' STORED AS TEXTFILE SELECT * FROM ctcf_repeats_cp;


-- Hi-C data

CREATE TABLE ctcf_repeats_hic ROW FORMAT DELIMITED FIELDS TERMINATED BY '\t' STORED AS RCFile AS SELECT cs.*, cp.id AS `id_cp`, cp.start1, cp.end1, cp.start2, cp.end2, cp.iab, cp.fdr FROM ctcf_repeats cs INNER JOIN chia_pet cp ON (cp.cell = cs.cell) WHERE cp.chrom1 = cs.chrom AND cp.chrom2 = cs.chrom AND cp.factor='Hi-C' AND ((cs.peak_ctcf BETWEEN cp.start1 AND cp.end1) OR (cs.peak_ctcf BETWEEN cp.start2 AND cp.end2));

CREATE TABLE ctcf_repeats_hic_mt ROW FORMAT DELIMITED FIELDS TERMINATED BY '\t' STORED AS RCFile AS SELECT cs.*, mt.id AS `id_mt`, mt.chromStart AS `chromStart_mt`, mt.chromEnd AS `chromEnd_mt`, mt.strand AS `strand_mt`, mt.score AS `score_mt` FROM ctcf_repeats_hic cs INNER JOIN motifs mt ON (mt.species = cs.species) WHERE cs.factor='CTCF' AND mt.factor = cs.factor AND mt.chrom = cs.chrom AND ( (mt.chromStart + mt.chromEnd)/2 BETWEEN cs.peak_ctcf-50 AND cs.peak_ctcf+50);

CREATE TABLE ctcf_hic ROW FORMAT DELIMITED FIELDS TERMINATED BY '\t' STORED AS RCFile AS SELECT cp.id AS `id_cp`, cp.chrom1 AS `chrom`, cp.start1, cp.end1, cp.start2, cp.end2, cp.iab, cp.fdr, cp.cell AS `cell`, cp.species AS `species`, cs.id AS `id_ctcf`, cs.chromstart AS `chromStart_ctcf`, cs.chromEnd AS `chromEnd_ctcf`, cs.signalValue AS `signalValue_ctcf`, cs.peak AS `peak_ctcf` FROM chia_pet cp INNER JOIN chip_seq cs ON (cs.cell = cp.cell) WHERE  cp.chrom1 = cs.chrom AND cp.chrom2 = cs.chrom AND cp.factor = 'Hi-C' AND cs.factor='CTCF' AND ((cs.peak BETWEEN cp.start1 AND cp.end1) OR (cs.peak BETWEEN cp.start2 AND cp.end2));

CREATE TABLE ctcf_hic_mt ROW FORMAT DELIMITED FIELDS TERMINATED BY '\t' STORED AS RCFile AS SELECT cp.*, mt.id AS `id_mt`, mt.chromstart AS `start_mt`, mt.chromend AS `end_mt`, mt.score AS `score_mt`, mt.strand AS `strand_mt` FROM ctcf_hic cp INNER JOIN motifs mt ON (mt.species = cp.species) WHERE mt.factor = 'CTCF' AND mt.chrom = cp.chrom AND ( (mt.chromStart + mt.chromEnd)/2 BETWEEN cp.peak_ctcf-50 AND cp.peak_ctcf+50);

INSERT OVERWRITE DIRECTORY '/user/adadiehl/ctcf_te_mt_intersection' ROW FORMAT DELIMITED FIELDS TERMINATED BY'\t' STORED AS TEXTFILE SELECT * FROM ctcf_repeats_hic_mt;
INSERT OVERWRITE DIRECTORY '/user/adadiehl/ctcf_mt_intersection' ROW FORMAT DELIMITED FIELDS TERMINATED BY'\t' STORED AS TEXTFILE SELECT * FROM ctcf_hic_mt;
INSERT OVERWRITE DIRECTORY '/user/adadiehl/ctcf_te_intersection' ROW FORMAT DELIMITED FIELDS TERMINATED BY'\t' STORED AS TEXTFILE SELECT * FROM ctcf_repeats_hic;
