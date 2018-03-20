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
	re_match = re.search('C(.+)?W[YLF][QL]', p_sequence_aa)
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

DROP FUNCTION IF EXISTS user_defined_functions.load_ab_aa_seq(TEXT, TEXT);
CREATE OR REPLACE FUNCTION user_defined_functions.load_ab_aa_seq(p_external_identifier TEXT, p_aa_seq TEXT)
RETURNS TABLE(aa_seq_md5 TEXT, external_identifier TEXT, aa_seq TEXT, aa_count INTEGER, seq_mol_wt REAL, cysteine_count INTEGER, 
			chain_type TEXT, cdr_1 TEXT, cdr_2 TEXT, cdr_3 TEXT)
AS
$$
DECLARE
  l_aa_seq_md5 TEXT := MD5(p_aa_seq);
  l_aa_seq TEXT := UPPER(TRIM(p_aa_seq));
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

SELECT aa_seq_md5, external_identifier, aa_seq, aa_count, seq_mol_wt, cysteine_count, chain_type, cdr_1, cdr_2, cdr_3
FROM
  user_defined_functions.load_ab_aa_seq('AB_H1', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGSINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARDGLMDVWGQGTAVTVSS');



