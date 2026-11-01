"""Temporary foreground progress monitor for this active analysis turn."""
from pathlib import Path
import json,time
P=Path(__file__).resolve().parent
def main():
    while True:
        try:
            q=json.loads((P/'qa/extraction_progress.json').read_text());lines=(P/'revision_pipeline.log').read_text(errors='replace').splitlines();stage=lines[-1] if lines else 'Waiting for complete cohort';print(time.strftime('%H:%M:%S'),'native fields',q['completed'],'/',q['planned'],'technical errors',q['technical_errors'],'|',stage,flush=True)
            if any('Initial version3 documents ready' in x or 'Revision step failed' in x or 'Manuscript build failed' in x for x in lines):break
        except (FileNotFoundError,json.JSONDecodeError):print('Progress checkpoint is being updated',flush=True)
        time.sleep(60)
if __name__=='__main__':main()
