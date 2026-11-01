from pathlib import Path
import pypdfium2 as pdfium
P=Path(__file__).resolve().parent/'qa'
pdf=pdfium.PdfDocument(str(P/'DYNLL1_U2OS_public_validation_sources.pdf'))
for i in range(len(pdf)):
    page=pdf[i];bitmap=page.render(scale=1.8);bitmap.to_pil().save(P/f'page-{i+1}.png');bitmap.close();page.close()
print('Rendered',len(pdf),'pages')
