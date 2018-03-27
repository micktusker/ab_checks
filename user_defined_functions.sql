CREATE OR REPLACE FUNCTION user_defined_functions.count_substrings(p_target_string TEXT, p_substring TEXT)
RETURNS INTEGER
AS
$$
	return p_target_string.count(p_substring)
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;

CREATE OR REPLACE FUNCTION user_defined_functions.first_cysteine_position(p_sequence_aa TEXT)
RETURNS INTEGER
AS
$$
	return p_sequence_aa.index('C') + 1
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;
SELECT user_defined_functions.first_cysteine_position('ABCDABCD');
CREATE OR REPLACE FUNCTION user_defined_functions.last_cysteine_position(p_sequence_aa TEXT)
RETURNS INTEGER
AS
$$
	return p_sequence_aa.rindex('C') + 1
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;
SELECT user_defined_functions.last_cysteine_position('ABCDABCD');


-- 'WYQ', 'WLQ', 'WFQ', 'WYL'
CREATE OR REPLACE FUNCTION user_defined_functions.get_cdr_l1_sequence(p_sequence_aa TEXT)
RETURNS TEXT
AS
$$
	import re
	re_match = re.search('C(.{10,})?W.[QL]', p_sequence_aa)
	if re_match:
		return re_match.groups()[0]
	else:
		return None
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;

--Start	always 16 residues after the end of L1
--Residues before 	generally Ile-Tyr, but also, Val-Tyr, Ile-Lys, Ile-Phe
--Length	always 7 residues (except NEW (7FAB) which has a deletion in this region)
CREATE OR REPLACE FUNCTION user_defined_functions.get_cdr_l2_sequence(p_sequence_aa TEXT, p_cdr_l1_sequence TEXT)
RETURNS TEXT
AS
$$
	import re
	if p_cdr_l1_sequence == None: return None
	pattern = p_cdr_l1_sequence + '.{15}(.{7})'
	re_match = re.search(pattern, p_sequence_aa)
	if re_match:
		return re_match.groups()[0]
	else:
		return None
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;

--Start	always 33 residues after end of L2 (except NEW (7FAB) which has the deletion at the end of CDR-L2)
--Residue before 	always Cys
--Residues after	always Phe-Gly-XXX-Gly
--Length	7 to 11 residues
--r'C.+C(.+)FG.G
CREATE OR REPLACE FUNCTION user_defined_functions.get_cdr_l3_sequence(p_sequence_aa TEXT)
RETURNS TEXT
AS
$$
	import re
	pattern = 'C.+C(.+)FG.G'
	re_match = re.search(pattern, p_sequence_aa)
	if re_match:
		return re_match.groups()[0]
	else:
		return None
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;

--CDR-H1
--Start	Approx residue 26 (always 4 after a Cys) [Chothia / AbM defintion];
--Kabat definition starts 5 residues later
--Residues before 	always Cys-XXX-XXX-XXX
--Residues after	always a Trp. Typically Trp-Val, but also, Trp-Ile, Trp-Ala
--Length	10 to 12 residues [AbM definition];
--Chothia definition excludes the last 4 residues
CREATE OR REPLACE FUNCTION user_defined_functions.get_cdr_h1_sequence(p_sequence_aa TEXT)
RETURNS TEXT
AS
$$
	import re
	pattern = 'C.{4}(.+?)W[VIA]'
	re_match = re.search(pattern, p_sequence_aa)
	if re_match:
		return re_match.groups()[0]
	else:
		return None
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;
--CDR-H2
--Start	always 15 residues after the end of Kabat / AbM definition) of CDR-H1
--Residues before 	typically Leu-Glu-Trp-Ile-Gly, but a number of variations
--Residues after	Lys/Arg-Leu/Ile/Val/Phe/Thr/Ala-Thr/Ser/Ile/Ala
--Length	Kabat definition 16 to 19 residues;
--AbM (and recent Chothia) definition ends 7 residues earlier
CREATE OR REPLACE FUNCTION user_defined_functions.get_cdr_h2_sequence(p_sequence_aa TEXT, p_cdr_h1_sequence TEXT, p_cdr_h3_sequence TEXT)
RETURNS TEXT
AS
$$
	import re
	if p_cdr_h1_sequence == None or p_cdr_h3_sequence == None: return None
	pattern = p_cdr_h1_sequence + '.{14}(.+).{32}' + p_cdr_h3_sequence
	re_match = re.search(pattern, p_sequence_aa)
	if re_match:
		return re_match.groups()[0]
	else:
		return None
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;
--CDR-H3
--Start	always 33 residues after end of CDR-H2 (always 2 after a Cys)
--Residues before 	always Cys-XXX-XXX (typically Cys-Ala-Arg)
--Residues after	always Trp-Gly-XXX-Gly
--Length	3 to 25(!) residues
CREATE OR REPLACE FUNCTION user_defined_functions.get_cdr_h3_sequence(p_sequence_aa TEXT)
RETURNS TEXT
AS
$$
	import re
	pattern = 'C.+C..(.+)WG.G'
	re_match = re.search(pattern, p_sequence_aa)
	if re_match:
		return re_match.groups()[0]
	else:
		return None
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;

CREATE OR REPLACE FUNCTION user_defined_functions.get_sequence_molecular_weight(p_sequence_aa TEXT)
RETURNS REAL
AS
$$
	from Bio.SeqUtils.ProtParam import ProteinAnalysis
	water_mol_wt = 18
	sequence_obj = ProteinAnalysis(p_sequence_aa)
	return round(sequence_obj.molecular_weight() + water_mol_wt, 2)
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;

CREATE OR REPLACE FUNCTION user_defined_functions.get_liability_doublet_info(p_sequence_aa TEXT)
RETURNS TEXT[]
AS
$$
	if p_sequence_aa == None: return None
	liability_doublets_found = []
	liability_doublets = ['NG', 'NM', 'NS', 'NT', 'DG', 'DS', 'DT', 'DD', 'DM']
	for liability_doublet in liability_doublets:
		if p_sequence_aa.find(liability_doublet) > -1:
			liability_doublets_found.append(liability_doublet)
	return liability_doublets_found
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;

DROP FUNCTION IF EXISTS user_defined_functions.load_ab_aa_seq(TEXT, TEXT, TEXT);
CREATE OR REPLACE FUNCTION user_defined_functions.load_ab_aa_seq(p_external_identifier TEXT, p_aa_seq TEXT, p_ab_target TEXT)
RETURNS TABLE(aa_seq_md5 TEXT, external_identifier TEXT, aa_seq TEXT, ab_target TEXT, aa_count INTEGER, seq_mol_wt REAL, cysteine_count INTEGER, 
			chain_type TEXT, cdr_1 TEXT, cdr_2 TEXT, cdr_3 TEXT)
AS
$$
DECLARE
  l_aa_seq_md5 TEXT := MD5(p_aa_seq);
  l_aa_seq TEXT := user_defined_functions.get_valid_amino_acid_sequence(p_aa_seq);
  l_ab_target TEXT := UPPER(TRIM(p_ab_target));
  l_aa_count INTEGER := LENGTH(l_aa_seq);
  l_seq_mol_wt REAL := user_defined_functions.get_sequence_molecular_weight(l_aa_seq);
  l_cysteine_count INTEGER := user_defined_functions.count_substrings(l_aa_seq, 'C');
  l_chain_type TEXT;
  l_cdr1 TEXT;
  l_cdr2 TEXT;
  l_cdr3 TEXT;
BEGIN
  IF l_aa_seq ~ 'TVSS$' THEN
    l_chain_type := 'H';
	l_cdr1 := user_defined_functions.get_cdr_h1_sequence(l_aa_seq);
	l_cdr3 := user_defined_functions.get_cdr_h3_sequence(l_aa_seq);
	l_cdr2 := user_defined_functions.get_cdr_h2_sequence(l_aa_seq, l_cdr1, l_cdr3);
  ELSIF l_aa_seq ~ 'EIK$' THEN
    l_chain_type := 'L';
	l_cdr1 := user_defined_functions.get_cdr_l1_sequence(l_aa_seq);
	l_cdr2 := user_defined_functions.get_cdr_l2_sequence(l_aa_seq, l_cdr1);
	l_cdr3 := user_defined_functions.get_cdr_l3_sequence(l_aa_seq);
  END IF;
  IF NOT EXISTS(SELECT 1 FROM ab_checker.antibody_sequences ab_seqs WHERE ab_seqs.aa_seq_md5 = l_aa_seq_md5) THEN
    INSERT INTO ab_checker.antibody_sequences(aa_seq_md5, 
											external_identifier, 
											aa_seq,
											ab_target,
											aa_count,
											seq_mol_wt, 
											cysteine_count, 
											chain_type, 
											cdr_1, 
											cdr_2, 
											cdr_3)
											VALUES(l_aa_seq_md5, 
												p_external_identifier, 
												l_aa_seq,
												l_ab_target,
												l_aa_count, 
												l_seq_mol_wt, 
												l_cysteine_count, 
												l_chain_type, 
												l_cdr1, 
												l_cdr2, 
												l_cdr3);
  END IF;
  RETURN QUERY
  SELECT 
    l_aa_seq_md5,
	p_external_identifier,
	l_aa_seq,
	l_ab_target,
	l_aa_count, 
	l_seq_mol_wt, 
	l_cysteine_count, 
	l_chain_type,
	l_cdr1,
	l_cdr2,
	l_cdr3;
END;
$$
LANGUAGE plpgsql
VOLATILE 
SECURITY DEFINER;

SELECT *
FROM
  user_defined_functions.load_ab_aa_seq('AB_H1', 'DIVMTQSPDSLAVSLGERATMSCKSSQSLLYSSNQKNYLAWHQQKPGQPPKLLIYWASTRESGVPDRFSGS GSGTDFTLTISSLQAEDLAIYYCQQYYTYPLTFGAGTKLEIK ', 'BAG3');

CREATE OR REPLACE FUNCTION user_defined_functions.get_valid_amino_acid_sequence(p_aa_seq TEXT)
RETURNS TEXT
AS
$$
	from Bio.Data import IUPACData
	import re
	if p_aa_seq == None: return None
	aa_seq = p_aa_seq.upper()
	aa_seq = re.sub(r'\s+', '', aa_seq)
	allowed_aa_set = set(list(IUPACData.protein_letters))
	given_sequence_aa_set = set(list(aa_seq))
	non_aa_set = given_sequence_aa_set - allowed_aa_set
	if len(non_aa_set) > 0: raise Exception('Given sequence contains non amino acid characters!')
	return aa_seq
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;
SELECT * FROM user_defined_functions.get_valid_amino_acid_sequence(' ac gt ');

DROP FUNCTION user_defined_functions.get_potential_aa_liabilty_info(TEXT);
CREATE OR REPLACE FUNCTION user_defined_functions.get_potential_aa_liabilty_info(p_aa_seq TEXT)
RETURNS TEXT[]
AS
$$
	aa_seq = str(plpy.execute("SELECT user_defined_functions.get_valid_amino_acid_sequence(%r)" % p_aa_seq))
	liabilities_found = []
	liabilities = ['NG', 'NM', 'NS', 'NT', 'DG', 'DS', 'DT', 'DD', 'DM', 'M', 'C']
	for liability in liabilities:
		if aa_seq.count(liability) > 0:
			liabilities_found.append(liability + ': ' + str(aa_seq.count(liability)))
	return liabilities_found
$$
LANGUAGE 'plpythonu'
STABLE
SECURITY DEFINER;

SELECT *
FROM
  user_defined_functions.get_potential_aa_liabilty_info('DIVMTQSPDSLAVSLGERATMSCKSSQSLLYSSNQKNYLAWHQQKPGQPPKLLIYWASTRESGVPDRFSGS GSGTDFTLTISSLQAEDLAIYYCQQYYTYPLTFGAGTKLEIK');

CREATE OR REPLACE FUNCTION user_defined_functions.get_cdr_liabilities(p_uniq_ab_id TEXT, p_cdr_name TEXT, p_cdr_seq TEXT)
RETURNS JSONB
AS
$$
DECLARE
  l_summary TEXT;
  l_liabilities TEXT[] := user_defined_functions.get_potential_aa_liabilty_info(p_cdr_seq);
BEGIN
  l_summary := FORMAT('{"cdr_name": "%s", "uniq_ab_id": "%s", "liabilities": [%s]}', 
					p_cdr_name, p_uniq_ab_id, '"' || 
					(ARRAY_TO_STRING(l_liabilities, '","')) || '"');
  RETURN l_summary::JSONB;
END
$$
LANGUAGE plpgsql
STABLE
SECURITY DEFINER;
SELECT * FROM user_defined_functions.get_cdr_liabilities('ab_xyz', 'cdr_l1', 'GMTFCSSLACIF');

CREATE OR REPLACE FUNCTION user_defined_functions.make_full_ab_id(p_full_l_chain_aa TEXT, p_full_h_chain_aa TEXT)
RETURNS TEXT
AS
$$
DECLARE
  l_full_l_chain_aa TEXT := REGEXP_REPLACE(UPPER(p_full_l_chain_aa), '\s', '', 'g');
  l_full_h_chain_aa TEXT := REGEXP_REPLACE(UPPER(p_full_h_chain_aa), '\s', '', 'g');
BEGIN
  RETURN MD5(l_full_l_chain_aa || l_full_h_chain_aa);
END;
$$
LANGUAGE plpgsql
IMMUTABLE
SECURITY DEFINER
RETURNS NULL ON NULL INPUT;

CREATE OR REPLACE FUNCTION user_defined_functions.get_pattern_location_info(p_text_to_search TEXT, p_pattern TEXT)
RETURNS JSONB
AS
$$
	import re
	import json
	pattern_location_matches = []
	match_info = {}
	for m in re.finditer(p_pattern, p_text_to_search):
		match_info['start'] = m.start()
		match_info['matched_text'] = m.group()
		match_info['end'] = m.end()
		pattern_location_matches.append(match_info)
	return json.dumps(pattern_location_matches)
$$
LANGUAGE plpythonu
IMMUTABLE
SECURITY DEFINER
RETURNS NULL ON NULL INPUT;
SELECT * FROM user_defined_functions.get_pattern_location_info('MNFGLRLIFLVLTLKGVQCQVQLVQSGAEVKKPGA', 'V.T');

