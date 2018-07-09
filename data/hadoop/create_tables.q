CREATE DATABASE IF NOT EXISTS ctcf_te_analysis_2;
CREATE EXTERNAL TABLE chia_pet(id INT, chrom1 STRING, start1 INT, end1 INT, chrom2 STRING,
                               start2 INT, end2 INT, iab INT, fdr FLOAT)
    PARTITIONED BY (factor STRING, species STRING, cell STRING)
    CLUSTERED BY (chrom1, chrom2) SORTED BY (start1, end1, start2, end2) INTO 32 BUCKETS
    ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'
    STORED AS TEXTFILE
    LOCATION '/user/adadiehl/db_store_2/chia_pet';

CREATE EXTERNAL TABLE chip_seq(id INT, chrom STRING, chromStart INT, chromEnd INT,
                               name STRING, score STRING, strand STRING, signalValue FLOAT,
                               pValue FLOAT, fdr FLOAT, peak INT)
    PARTITIONED BY (factor STRING, species STRING, cell STRING)
    CLUSTERED BY (chrom) SORTED BY (chromStart, chromEnd) INTO 32 BUCKETS
    ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'
    STORED AS TEXTFILE
    LOCATION '/user/adadiehl/db_store_2/chip_seq';

CREATE EXTERNAL TABLE repeatmasker(id INT, chrom STRING, chromStart INT, chromEnd INT,
                                   name STRING, class STRING, family STRING, pctDiv FLOAT)
    PARTITIONED BY (species STRING)
    CLUSTERED BY (chrom) SORTED BY (chromStart, chromEnd) INTO 32 BUCKETS
    ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'
    STORED AS TEXTFILE
    LOCATION '/user/adadiehl/db_store_2/repeatmasker';

CREATE EXTERNAL TABLE motifs(id INT, chrom STRING, chromStart INT, chromEnd INT,
                             strand STRING, score FLOAT)
    PARTITIONED BY (factor STRING, species STRING)
    CLUSTERED BY (chrom) SORTED BY (chromStart, chromEnd) INTO 32 BUCKETS
    ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'
    STORED AS TEXTFILE
    LOCATION '/user/adadiehl/db_store_2/motifs';

CREATE EXTERNAL TABLE chia_pet_st(id INT, chrom1 STRING, start1 INT, end1 INT, chrom2 STRING,
                                  start2 INT, end2 INT, iab INT, fdr FLOAT, factor STRING,
				  species STRING, cell STRING)
    ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'
    STORED AS TEXTFILE
    LOCATION '/user/adadiehl/db_store_2/chia_pet';

CREATE EXTERNAL TABLE chip_seq_st(id INT, chrom STRING, chromStart INT, chromEnd INT,
                                  name STRING, score STRING, strand STRING, signalValue FLOAT,
                                  pValue FLOAT, fdr FLOAT, peak INT, factor STRING,
				  species STRING, cell STRING)
    ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'
    STORED AS TEXTFILE
    LOCATION '/user/adadiehl/db_store_2/chip_seq';

CREATE EXTERNAL TABLE repeatmasker_st(id INT, chrom STRING, chromStart INT, chromEnd INT,
                                      name STRING, class STRING, family STRING, pctDiv FLOAT,
				      species STRING)
    ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'
    STORED AS TEXTFILE
    LOCATION '/user/adadiehl/db_store_2/repeatmasker';

CREATE EXTERNAL TABLE motifs_st(id INT, chrom STRING, chromStart INT, chromEnd INT,
                                strand STRING, score FLOAT, factor STRING, species STRING)
    ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'
    STORED AS TEXTFILE
    LOCATION '/user/adadiehl/db_store_2/motifs';

LOAD DATA INPATH '/user/adadiehl/chia_pet.dat' OVERWRITE INTO TABLE chia_pet_st;

INSERT INTO chia_pet PARTITION (factor='RAD21', species='hg19', cell='GM12878') SELECT id, chrom1, start1, end1, chrom2, start2, end2, iab, fdr FROM chia_pet_st WHERE factor='RAD21' AND species='hg19' AND cell='GM12878';
INSERT INTO chia_pet PARTITION (factor='RAD21', species='hg19', cell='K562') SELECT id, chrom1, start1, end1, chrom2, start2, end2, iab, fdr FROM chia_pet_st WHERE factor='RAD21' AND species='hg19' AND cell='K562';
INSERT INTO chia_pet PARTITION (factor='Hi-C', species='hg19', cell='GM12878') SELECT id, chrom1, start1, end1, chrom2, start2, end2, iab, fdr FROM chia_pet_st WHERE factor='Hi-C' AND species='hg19' AND cell='GM12878';
INSERT INTO chia_pet PARTITION (factor='Hi-C', species='hg19', cell='K562') SELECT id, chrom1, start1, end1, chrom2, start2, end2, iab, fdr FROM chia_pet_st WHERE factor='Hi-C' AND species='hg19' AND cell='K562';
INSERT INTO chia_pet PARTITION (factor='Hi-C', species='mm9', cell='CH12') SELECT id, chrom1, start1, end1, chrom2, start2, end2, iab, fdr FROM chia_pet_st WHERE factor='Hi-C' AND species='mm9' AND cell='CH12';


LOAD DATA INPATH '/user/adadiehl/chip_seq.dat' OVERWRITE INTO TABLE chip_seq_st;
INSERT INTO chip_seq PARTITION (factor='CTCF', species='hg19', cell='GM12878') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='CTCF' AND species='hg19' AND cell='GM12878';
INSERT INTO chip_seq PARTITION (factor='RAD21', species='hg19', cell='GM12878') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='RAD21' AND species='hg19' AND cell='GM12878';
INSERT INTO chip_seq PARTITION (factor='SMC3', species='hg19', cell='GM12878') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='SMC3' AND species='hg19' AND cell='GM12878';

INSERT INTO chip_seq PARTITION (factor='CTCF', species='hg19', cell='K562') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='CTCF' AND species='hg19' AND cell='K562';
INSERT INTO chip_seq PARTITION (factor='RAD21', species='hg19', cell='K562') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='RAD21' AND species='hg19' AND cell='K562';
INSERT INTO chip_seq PARTITION (factor='SMC3', species='hg19', cell='K562') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='SMC3' AND species='hg19' AND cell='K562';

INSERT INTO chip_seq PARTITION (factor='CTCF', species='mm9', cell='CH12') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='CTCF' AND species='mm9' AND cell='CH12';
INSERT INTO chip_seq PARTITION (factor='RAD21', species='mm9', cell='CH12') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='RAD21' AND species='mm9' AND cell='CH12';
INSERT INTO chip_seq PARTITION (factor='SMC3', species='mm9', cell='CH12') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='SMC3' AND species='mm9' AND cell='CH12';

INSERT INTO chip_seq PARTITION (factor='CTCF', species='mm9', cell='MEL') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='CTCF' AND species='mm9' AND cell='MEL';
INSERT INTO chip_seq PARTITION (factor='RAD21', species='mm9', cell='MEL') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='RAD21' AND species='mm9' AND cell='MEL';
INSERT INTO chip_seq PARTITION (factor='SMC3', species='mm9', cell='MEL') SELECT id, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, fdr, peak FROM chip_seq_st WHERE factor='SMC3' AND species='mm9' AND cell='MEL';


LOAD DATA INPATH '/user/adadiehl/repeatmasker.dat' OVERWRITE INTO TABLE repeatmasker_st;
INSERT INTO repeatmasker PARTITION (species='hg19') SELECT id, chrom, chromStart, chromEnd, name, class, family, pctDiv FROM repeatmasker_st WHERE species='hg19';
INSERT INTO repeatmasker PARTITION (species='mm9') SELECT id, chrom, chromStart, chromEnd, name, class, family, pctDiv FROM repeatmasker_st WHERE species='mm9';


LOAD DATA INPATH '/user/adadiehl/motifs.dat' OVERWRITE INTO TABLE motifs_st;
INSERT INTO motifs PARTITION (species='hg19', factor='CTCF') SELECT id, chrom, chromStart, chromEnd, strand, score FROM motifs_st WHERE species='hg19' AND factor='CTCF';
INSERT INTO motifs PARTITION (species='mm9', factor='CTCF') SELECT id, chrom, chromStart, chromEnd, strand, score FROM motifs_st WHERE species='mm9' AND factor='CTCF';
