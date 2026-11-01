from pathlib import Path
import pypdfium2 as pdfium
P=Path(__file__).resolve().parent/'qa/render'
pdf=pdfium.PdfDocument(str(P/'manuscript.pdf'))
for i in range(len(pdf)):
    page=pdf[i];bitmap=page.render(scale=1.6);bitmap.to_pil().save(P/f'page-{i+1}.png');bitmap.close();page.close()
print('Rendered',len(pdf),'pages')
