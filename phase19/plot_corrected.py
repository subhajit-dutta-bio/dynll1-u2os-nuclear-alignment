"""Redraw the three corrected figures and update their existing source registry."""
from figure_common import *
from plot_supp import supp2,supp3,supp5
def main():
    prior=json.loads((ROOT/'phase17/panel_registry.json').read_text())
    d=Data()
    for fn in [supp2,supp3,supp5]:fn(d)
    changed={(r['figure'],r['panel']):r for r in REG};assert len(changed)==24
    current=[]
    for r in prior:
        r=changed.get((r['figure'],r['panel']),r.copy())
        r['source_data']=r['source_data'].replace('phase17/panel_data','phase19/panel_data');current.append(r)
    assert len(current)==88 and len({(r['figure'],r['panel']) for r in current})==88
    (OUT/'panel_registry.json').write_text(json.dumps(current,indent=2),encoding='utf-8');pd.DataFrame(current).to_csv(OUT/'panel_registry.csv',index=False)
    print('Redrew Supplementary Figures 2, 3 and 5; preserved all 88 panel records.')
if __name__=='__main__':main()
