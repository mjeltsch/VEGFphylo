BEGIN TRANSACTION;
CREATE TABLE IF NOT EXISTS `species` (
	`scientific_name`	TEXT NOT NULL,
	`taxon_id`	INTEGER NOT NULL,
	`phylum`	TEXT NOT NULL,
	PRIMARY KEY(`taxon_id`)
);
CREATE TABLE IF NOT EXISTS `protein` (
	`gid`	INTEGER NOT NULL,
	`protein_id`	TEXT,
	`fasta_description`	TEXT NOT NULL,
	`species`	TEXT NOT NULL,
	`ortholog_group`	TEXT,
	`curated_manually_by`	TEXT DEFAULT ('NULL'),
	PRIMARY KEY(`gid`)
);
CREATE TABLE IF NOT EXISTS `ortholog_groups` (
	`ortholog_group`	TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS `curator` (
	`curator`	TEXT NOT NULL
);
COMMIT;
