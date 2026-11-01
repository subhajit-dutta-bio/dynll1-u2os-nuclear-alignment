"""Bibliography corrected against primary publication and resource records."""
from pathlib import Path
import json
P=Path(__file__).resolve().parent;ROOT=P.parent
def references():
    r={int(k):v for k,v in json.loads((ROOT/'phase10/data/crossref_references.json').read_text(encoding='utf-8')).items()}
    # Repair missing sentence boundaries in the retained registry, without
    # changing source titles or treating the registry as raw Crossref JSON.
    for i,name in [(2,'Journal of Biological Chemistry'),(3,'Nucleic Acids Research'),(4,'Nature.'),(5,'Nature Structural'),(7,'Nature Cell Biology'),(8,'Nature Chemical Biology'),(9,'Nature Communications'),(13,'Nature.'),(15,'The EMBO Reports'),(16,'Nature Communications')]:
        r[i]=r[i].replace(' '+name,'. '+name,1)
    r[15]=r[15].replace('The EMBO Reports','EMBO Reports').replace('2019;20(9).','2019;20(9):e47495.')
    r.update({
    1:'1. Chandrasekaran SN et al. Morphological map of under- and overexpression of genes in human cells. Nature Methods. 2025;22:1742–1752. https://doi.org/10.1038/s41592-025-02753-9. Image collection cpg0016-jump.',
    6:'6. Broad Bioimage Benchmark Collection. BBBC039v1 annotated U2OS nuclear images and official reference splits. https://bbbc.broadinstitute.org/BBBC039. Accessed 30 October 2026.',
    9:'9. Allen TP, Roennfeldt AE, Reckdharajkumar M, Sullivan AE, Liu M, Quinn RJ et al. dFLASH; dual FLuorescent transcription factor activity sensor for histone integrated live-cell reporting and high-content screening. Nature Communications. 2025;16:3298. https://doi.org/10.1038/s41467-025-58488-w. Counts GSE290538.',
    10:'10. Harper NW, Birdsall GA, Honeywell ME, Ward KM, Pai AA, Lee MJ. RNA Pol II inhibition activates cell death independently from the loss of transcription. Cell. 2025;188(22):6301–6316.e29. https://doi.org/10.1016/j.cell.2025.07.034. Counts GSE283147.',
    11:'11. McDougal MB, DeMaria AM, Ohlson MB, Kumar A, Xing C, Schoggins JW. Interferon inhibits a model RNA virus via a limited set of inducible effector genes. EMBO Reports. 2023;24(9):e56901. https://doi.org/10.15252/embr.202356901. Counts GSE233195.',
    12:'12. DepMap Consortium. DepMap 24Q4 Public. https://doi.org/10.25452/figshare.plus.27993248. U2OS model ACH-000364 and screen SC-000364.AV01.',
    14:'14. Cell Maps for Artificial Intelligence. U2OS protein-interaction release. https://doi.org/10.18130/V3/164SI1. Associated MassIVE resource MSV000097168 and NDEx network 95bc75d5-d1d1-11ee-8a40-005056ae23aa.',
    17:'17. Caicedo JC et al. Evaluation of Deep Learning Strategies for Nucleus Segmentation in Fluorescence Images. Cytometry Part A. 2019;95(9):952–965. https://doi.org/10.1002/cyto.a.23863.',
    18:'18. HUGO Gene Nomenclature Committee. HGNC complete gene set. https://www.genenames.org/download/archive/. Retained complete-set annotation retrieved 30 October 2026.',
    19:'19. Pachitariu M, Rariden M, Stringer C. Cellpose-SAM: superhuman generalization for cellular segmentation. bioRxiv. 2025. https://doi.org/10.1101/2025.04.28.651001. Implementation Cellpose4.2.1.1 with the retained cpsam_v2 checkpoint; https://github.com/MouseLand/cellpose.',
    20:'20. Lazar NH, Celik S, Chen L, Fay MM, Irish JC, Jensen J et al. High-resolution genome-wide mapping of chromosome-arm-scale truncations induced by CRISPR–Cas9 editing. Nature Genetics. 2024;56(7):1482–1493. https://doi.org/10.1038/s41588-024-01758-y.',
    21:'21. Okumura M, Natsume T, Kanemaki MT, Kiyomitsu T. Dynein–Dynactin–NuMA clusters generate cortical spindle-pulling forces as a multi-arm ensemble. eLife. 2018;7:e36559. https://doi.org/10.7554/eLife.36559.',
    22:'22. Saw TB, Doostmohammadi A, Nier V et al. Topological defects in epithelia govern cell death and extrusion. Nature. 2017;544:212–216. https://doi.org/10.1038/nature21718.',
    23:'23. Splinter D, Tanenbaum ME, Lindqvist A, Jaarsma D, Flotho A, Yu KL et al. Bicaudal D2, Dynein, and Kinesin-1 Associate with Nuclear Pore Complexes and Regulate Centrosome and Nuclear Positioning during Mitotic Entry. PLOS Biology. 2010;8(4):e1000350. https://doi.org/10.1371/journal.pbio.1000350.',
    24:'24. Wong CH, Wingett SW, Qian C, Hunter MR, Taliaferro JM, Ross-Thriepland D, Bullock SL. Genome-scale requirements for dynein-based transport revealed by a high-content arrayed CRISPR screen. Journal of Cell Biology. 2024;223(5):e202306048. https://doi.org/10.1083/jcb.202306048.',
    25:'25. Broad Genetic Perturbation Platform and Addgene. Human genome-wide CRISPR knockout Brunello library and official guide annotations. https://www.addgene.org/pooled-library/broadgpp-human-knockout-brunello/. Accessed 30 October 2026. Median intended cut positions are used only as gene-position proxies in the JUMP comparator annotation.'})
    r[16]=r[16].replace('2026;17(1).','2026;17:6857.')
    r[7]+=' Counts GSE226994. Author correction https://doi.org/10.1038/s41556-026-01884-1 concerns the Extended Data Figure8b image; no changed count table was identified.'
    r[8]+=' Counts GSE232592.'
    for i in [7,8]:r[i]=r[i].replace(' Counts GSE','. Counts GSE',1)
    return [r[i] for i in range(1,26)]
