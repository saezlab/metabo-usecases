# Load swisslipid
import importlib
import lipinet.parse_swisslipids
importlib.reload(lipinet.parse_swisslipids)  # reload the module after edits
from lipinet.parse_swisslipids import parse_swisslipids_data

sl_results = parse_swisslipids_data(verbose=False, use_cache=True)
df_sl_nodes = sl_results['df_nodes']
df_sl_edges = sl_results['df_edges']

# Check lipids by category
# CE
df_sl_nodes['Abbreviation_cleaned'] = df_sl_nodes['Abbreviation*'].str.replace(r'\([0-9]+[EZ](,[0-9]+[EZ])*\)', '', regex=True)
morita_ce = data[data['CompoundName'].str.startswith('CE(', na=False)]['CompoundName'].tolist()
for ce in morita_ce:
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == ce]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{ce}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# Cer
morita_cer = data[data['CompoundName'].str.startswith('Cer ', na=False)]['CompoundName']
morita_cer_converted = morita_cer.str.replace('Cer ', 'Cer(').str.replace('_', '/') + ')'
for orig, converted in zip(morita_cer, morita_cer_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# DG
morita_dg = data[data['CompoundName'].str.startswith('DG ', na=False)]['CompoundName']
morita_dg_converted = morita_dg.str.replace('DG ', 'DG(') + ')'
for orig, converted in zip(morita_dg, morita_dg_converted):
    rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
    chebi = rows['CHEBI'].dropna().unique()
    print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# FA
morita_fa = data[data['CompoundName'].str.startswith('FA(', na=False)]['CompoundName']
for fa in morita_fa:
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == fa]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{fa}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# Hexcer
morita_hexcer = data[data['CompoundName'].str.startswith('HexCer ', na=False)]['CompoundName']
morita_hexcer_converted = morita_hexcer.str.replace('HexCer ', 'HexCer(').str.replace('_', '/') + ')'
for orig, converted in zip(morita_hexcer, morita_hexcer_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# LPC
morita_lpc = data[data['CompoundName'].str.startswith('LPC ', na=False)]['CompoundName']
morita_lpc_converted = morita_lpc.str.replace('LPC ', 'LPC(') + ')'
for orig, converted in zip(morita_lpc, morita_lpc_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# LPE
morita_lpe = data[data['CompoundName'].str.startswith('LPE ', na=False)]['CompoundName']
morita_lpe_converted = morita_lpe.str.replace('LPE ', 'LPE(') + ')'
for orig, converted in zip(morita_lpe, morita_lpe_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# MG
morita_mg = data[data['CompoundName'].str.startswith('MG ', na=False)]['CompoundName']
morita_mg_converted = morita_mg.str.replace('MG ', 'MG(') + ')'

for orig, converted in zip(morita_mg, morita_mg_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# PA
morita_pa = data[data['CompoundName'].str.startswith('PA ', na=False)]['CompoundName']
morita_pa_converted = morita_pa.str.replace('PA ', 'PA(') + ')'

for orig, converted in zip(morita_pa, morita_pa_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# PC
morita_pc = data[data['CompoundName'].str.startswith('PC ', na=False)]['CompoundName']
morita_pc_converted = morita_pc.str.replace('PC ', 'PC(') + ')'

for orig, converted in zip(morita_pc, morita_pc_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# PE
morita_pe = data[data['CompoundName'].str.startswith('PE ', na=False)]['CompoundName']
morita_pe_converted = morita_pe.str.replace('PE ', 'PE(') + ')'
for orig, converted in zip(morita_pe, morita_pe_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# PG
morita_pg = data[data['CompoundName'].str.startswith('PG ', na=False)]['CompoundName']
morita_pg_converted = morita_pg.str.replace('PG ', 'PG(') + ')'

for orig, converted in zip(morita_pg, morita_pg_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# PI
morita_pi = data[data['CompoundName'].str.startswith('PI ', na=False)]['CompoundName']
morita_pi_converted = morita_pi.str.replace('PI ', 'PI(') + ')'

for orig, converted in zip(morita_pi, morita_pi_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# PS
morita_ps = data[data['CompoundName'].str.startswith('PS ', na=False)]['CompoundName']
morita_ps_converted = morita_ps.str.replace('PS ', 'PS(') + ')'

for orig, converted in zip(morita_ps, morita_ps_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# SM
morita_sm = data[data['CompoundName'].str.startswith('SM ', na=False)]['CompoundName']
morita_sm_converted = morita_sm.str.replace('SM ', 'SM(d') + ')'

for orig, converted in zip(morita_sm, morita_sm_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")

# TG
morita_tg = data[data['CompoundName'].str.startswith('TG(', na=False)]['CompoundName']
morita_tg_converted = morita_tg.str.replace(')(', '_', regex=False)

for orig, converted in zip(morita_tg, morita_tg_converted):
      rows = df_sl_nodes[df_sl_nodes['Abbreviation_cleaned'] == converted]
      chebi = rows['CHEBI'].dropna().unique()
      print(f"{orig}: {list(chebi) if len(chebi) > 0 else 'No CHEBI'}")
