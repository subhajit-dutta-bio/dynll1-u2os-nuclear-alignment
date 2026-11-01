from figure_common import *
from plot_main import main_figures
from plot_supp import supplementary_figures
def main():
    assert json.loads((P/'qa/numeric_verification.json').read_text())['all_pass'];d=Data();main_figures(d);supplementary_figures(d);assert len(REG)==80 and len({(x['figure'],x['panel']) for x in REG})==80
    (P/'panel_registry.json').write_text(json.dumps(REG,indent=2),encoding='utf-8');pd.DataFrame(REG).to_csv(P/'panel_registry.csv',index=False);print('Rendered ten figures and all80 numeric panel sources')
if __name__=='__main__':main()
