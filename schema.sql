DROP TABLE IF EXISTS ab_checker.antibody_sequences;
CREATE TABLE ab_checker.antibody_sequences(
  aa_seq_md5 TEXT PRIMARY KEY,
  external_identifier TEXT NOT NULL,
  aa_seq TEXT,
  ab_target TEXT,
  aa_count INTEGER,
  seq_mol_wt REAL,
  cysteine_count INTEGER,
  chain_type TEXT CHECK(chain_type IN('L', 'H')),
  cdr_1 TEXT,
  cdr_2 TEXT,
  cdr_3 TEXT
);


CREATE OR REPLACE VIEW ab_checker.vw_extracted_cdrs AS
SELECT
  sequence_id,
  sequence_aa,
  user_defined_functions.get_sequence_molecular_weight(sequence_aa) mol_wt,
  LENGTH(sequence_aa) amino_acid_count,
  CASE
    WHEN LENGTH(sequence_aa) <= 110 THEN 'L'
	ELSE 'H'
  END chain_type,
  sequence_md5,
  user_defined_functions.count_substrings(sequence_aa, 'C') cysteine_count,
  user_defined_functions.count_substrings(sequence_aa, 'M') methionine_count,
  user_defined_functions.first_cysteine_position(sequence_aa),
  user_defined_functions.last_cysteine_position(sequence_aa),
  user_defined_functions.get_cdr_l1_sequence(sequence_aa) cdr_l1,
  user_defined_functions.get_cdr_l2_sequence(sequence_aa, user_defined_functions.get_cdr_l1_sequence(sequence_aa)) cdr_l2,
  user_defined_functions.get_cdr_l3_sequence(sequence_aa) cdr_l3,
  user_defined_functions.get_cdr_h1_sequence(sequence_aa) cdr_h1,
  user_defined_functions.get_cdr_h2_sequence(sequence_aa, user_defined_functions.get_cdr_h1_sequence(sequence_aa), user_defined_functions.get_cdr_h3_sequence(sequence_aa)) cdr_h2,
  user_defined_functions.get_cdr_h3_sequence(sequence_aa) cdr_h3
FROM
  ab_checker.input_sequence
ORDER BY
  chain_type;

CREATE OR REPLACE VIEW ab_checker.vw_liability_doublet_info AS
SELECT
  sequence_md5,
  user_defined_functions.get_liability_doublet_info(cdr_l1) cdr_l1_liabilities,
  user_defined_functions.get_liability_doublet_info(cdr_l2) cdr_l2_liabilities,
  user_defined_functions.get_liability_doublet_info(cdr_l3) cdr_l3_liabilities,
  user_defined_functions.get_liability_doublet_info(cdr_h1) cdr_h1_liabilities,
  user_defined_functions.get_liability_doublet_info(cdr_h2) cdr_h2_liabilities,
  user_defined_functions.get_liability_doublet_info(cdr_h3) cdr_h3_liabilities
FROM
  ab_checker.vw_extracted_cdrs;
