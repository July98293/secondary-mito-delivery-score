

import os, sys, io, warnings
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

import xlrd
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, "data")
OUT  = os.path.join(ROOT, "otput")


MODULES = {
    "Exosome Machinery": {
        "weight": 0.20, "category": "core",
        "genes": {
            "TSG101":  "Tumor Susceptibility Gene 101; componente ESCRT-I, riconosce l'ubiquitina sui cargo degli MVB",
            "PDCD6IP": "ALIX; adattatore ESCRT-III, lega LBPA e orquestra la gemmazione di esosomi",
            "VPS4A":   "ATPasi che disassembla i complessi ESCRT-III sulla membrana dei MVB",
            "VPS4B":   "Isoforma di VPS4A; coopera con VPS4A nel rimodellamento delle membrane MVB",
            "CD63":    "Tetraspanica LAMP-3; marcatore canonico di esosomi, presente su lisosomi e MVB",
            "CD9":     "Tetraspanica; facilita la fusione di membrane e l'ancoraggio di esosomi alle cellule riceventi",
            "CD81":    "Tetraspanica; forma complessi con CD9/CD63, stabilizza i domini lipidici degli esosomi",
            "SDCBP":   "Sindecano-binding protein (sintenina-1); ponticella ESCRT e syndecan nell'endosorting",
        }
    },
    "Donor Fitness (Mito Health)": {
        "weight": 0.18, "category": "core",
        "genes": {
            "TFAM":     "Mitochondrial Transcription Factor A; compatta il mtDNA e ne regola la replicazione/trascrizione",
            "OPA1":     "GTPasi della fusione della membrana interna mitocondriale; mantiene le creste e il potenziale di membrana",
            "MFN1":     "Mitofusina-1; GTPasi della fusione della membrana esterna; coopera con MFN2",
            "MFN2":     "Mitofusina-2; fusione OMM e ancoraggio ER-mito; mutata in Charcot-Marie-Tooth 2A",
            "DNM1L":    "DRP1 (Dynamin-related protein 1); GTPasi della fissione mitocondriale",
            "PINK1":    "PTEN-Induced Kinase 1; sensore del potenziale di membrana; fosforila Parkin per mitofagia selettiva",
            "PRKN":     "Parkina (E3 ubiquitin ligase); ubiquitina proteine OMM danneggiate, segnale per mitofagia",
            "PPARGC1A": "PGC-1α; coattivatore trascrizionale principale della biogenesi mitocondriale",
            "ATP5F1A":  "Subunità α dell'ATP-sintasi (complesso V); catalizza la sintesi di ATP",
        }
    },
    "Rab Trafficking": {
        "weight": 0.12, "category": "secondary",
        "genes": {
            "RAB27A": "GTPasi Rab; regola l'aggancio dei granuli secretori/MVB alla membrana plasmatica (via melanofillina)",
            "RAB27B": "Isoforma di RAB27A; predominante in mastociti e basofili per la via secretoria",
            "RAB11A": "Regola il riciclo endosomale lento; coinvolto nel sorting di esosomi da compartimenti di riciclo",
            "RAB11B": "Isoforma di RAB11A; espresso prevalentemente nel cervello e nella muscolatura",
            "RAB35":  "Controllo del riciclo rapido endosomale e della secrezione di esosomi in modo ARF6-dipendente",
            "RAB5A":  "GTPasi degli endosomi precoci; regola la fusione omo- ed eterotipica tra endosomi",
            "RAB7A":  "Maturazione degli endosomi tardivi; trafficking verso lisosomi e biogenesi dei MVB",
        }
    },
    "MDV Pathway (Mitochondria-Derived Vesicles)": {
        "weight": 0.10, "category": "secondary",
        "genes": {
            "VPS35":  "Subunità cargo-recognition del Retromer; recupera proteine OMM nei MDV verso gli endosomi",
            "SNX9":   "Sorting Nexin 9; sente la curvatura di membrana e recluta la macchineria di fissione per i MDV",
            "ATG9A":  "Unica proteina transmembrana del pathway autofagico; veicola lipidi per la formazione di fagofori e MDV",
            "MUL1":   "Mitochondrial Ubiquitin Ligase 1 (MAPL); ubiquitina proteine OMM per il sorting nei MDV",
            "RAB9A":  "Trafficking dal compartimento tardo-endosomale al TGN; necessario per la maturazione dei MDV",
            "OPTN":   "Optineurina; adattatore ubiquitina-dipendente per mitofagia selettiva e MDV verso il Golgi",
        }
    },
    "TNT Capacity (Tunnelling Nanotubes)": {
        "weight": 0.08, "category": "secondary",
        "genes": {
            "RHOA":  "GTPasi Rho; promuove la polimerizzazione di actina e la retrazione di membrana necessaria per TNT",
            "CDC42": "GTPasi Cdc42; induce proiezioni di membrana filopodiali, precursori dei nanotubuli",
            "ACTA2": "α-actina del muscolo liscio; componente del citoscheletro contrattile che sostiene i TNT",
            "INF2":  "Formina che accelera la polimerizzazione e depolimerizzazione di actina; essenziale per TNT in cellule T",
            "FMNL1": "Formin-like 1; regola la polimerizzazione di actina alle proiezioni di membrana nelle cellule NK",
            "FMNL2": "Formin-like 2; promuove la lamellipodio/filipodio nelle cellule migratorie",
            "ARPC2": "Subunità del complesso Arp2/3; ramificazione del filamento di actina per la formazione di TNT",
            "ARPC3": "Subunità del complesso Arp2/3; coopera con ARPC2 nel nucleare nuovi filamenti",
        }
    },
    "Adhesion & Migration": {
        "weight": 0.08, "category": "secondary",
        "genes": {
            "CXCR4":  "Recettore per SDF-1/CXCL12; guida l'homing al midollo osseo e ai tessuti ipossici",
            "CCR7":   "Recettore per CCL19/21; dirige la migrazione ai linfonodi e alle zone T",
            "ITGA4":  "Integrina α4 (CD49d); si lega a VCAM-1 e fibronectina, media il rolling e l'adesione vascolare",
            "ITGB1":  "Integrina β1 (CD29); forma eterodimeri con >10 subunità α; adesione a ECM e ligandi cellulari",
            "ICAM1":  "ICAM-1 (CD54); ligando di LFA-1; stabilizza la sinapsi immunologica tra cellule donatrici e riceventi",
        }
    },
    "Metabolic Vigor": {
        "weight": 0.08, "category": "secondary",
        "genes": {
            "PFKFB3": "Bisfosfofruttocinasi; produce fruttosio-2,6-bisfosfato, attivatore allosterico della glicolisi",
            "SLC2A1": "GLUT1; trasportatore di glucosio basolaterale ubiquitario; aumenta la captazione energetica",
            "LDHA":   "Lattato deidrogenasi A; converte piruvato in lattato; sostiene la glicolisi anaerobica",
            "CPT1A":  "Carnitine Palmitoyltransferase 1A; enzima limitante per l'importazione degli acidi grassi nel mitocondrio",
            "HIF1A":  "HIF-1α; fattore di trascrizione dell'ipossia; up-regola GLUT1, LDHA e VEGF",
        }
    },
    "Microvesicle Formation": {
        "weight": 0.08, "category": "secondary",
        "genes": {
            "ARF6":   "GTPasi ARF6; regola il riciclaggio endosomale e la gemmazione di microvescicole dalla membrana plasmatica",
            "ROCK1":  "Rho-kinase 1; fosforila la catena leggera della miosina, contrae l'actina-cortex per la gemmazione",
            "ROCK2":  "Rho-kinase 2; isoforma di ROCK1 con funzioni simili nella contrazione del citoscheletro",
            "MYH9":   "Catena pesante della miosina IIA non-muscolare; motore molecolare per la gemmazione di MV",
            "ARRDC1": "Arrestin Domain-Containing 1; recluta TSG101 alla membrana plasmatica per gemmazione di MV TSG101+",
            "FLOT1":  "Flotillina-1; componente dei lipid raft; scaffolding per la gemmazione di MV di membrana plasmatica",
            "FLOT2":  "Flotillina-2; forma eterotetrameri con FLOT1, necessaria per organizzare i raft di membrana",
            "PLD2":   "Fosfolipasi D2; produce acido fosfatidico, promuove la curvatura di membrana per la gemmazione",
        }
    },
    "Immune Shedding (Safety+)": {
        "weight": 0.04, "category": "safety",
        "genes": {
            "CD63":  "→ vedi Exosome Machinery; anche marcatore di degranulazione controllata e sicurezza di secrezione",
            "LAMP1": "Proteina di membrana lisosomale 1 (CD107a); esportata sulla superficie durante la degranulazione; marker di secrezione controllata",
            "TSG101":"→ vedi Exosome Machinery; la sua presenza in MV indica un pathway ESCRT-dipendente (sicuro)",
        }
    },
    "Activation Markers (Safety−)": {
        "weight": 0.04, "category": "safety",
        "genes": {
            "IL2RA": "CD25; catena α del recettore per IL-2; upregolata in cellule T attivate e Treg; rischio infiammatorio",
            "CD69":  "Primo marcatore di attivazione cellulare; trattiene i linfociti nei tessuti (downregola S1PR1)",
            "TNF":   "TNF-α; citochina pro-infiammatoria sistemica; rischio principale nella somministrazione terapeutica",
            "IL1B":  "IL-1β; pirogeno endogeno e attivatore dell'inflammasoma; rischio di risposta infiammatoria sistemica",
            "IFNG":  "IFN-γ; citochina Th1; attiva macrofagi e upregola HLA-II; rischio immunogenicità",
            "HLA-DRA": "HLA-DR α; presentazione antigenica MHC classe II; indica cellule APC con potenziale immunogenico",
        }
    },
    "Viability / Survival": {
        "weight": None, "category": "viability",
        "genes": {
            "BCL2":  "BCL-2; proto-oncogene anti-apoptotico; stabilizza la membrana mitocondriale esterna, inibisce BAX",
            "MCL1":  "Myeloid Cell Leukemia 1; anti-apoptotico a breve emivita; critico per la sopravvivenza dei granulociti",
            "BIRC2": "Inhibitor of Apoptosis 2 (cIAP1); ubiquitina RIPK1 e blocca la via estrinseca dell'apoptosi",
            "XIAP":  "X-linked IAP; inibitore diretto delle caspasi 3, 7 e 9; gene di sopravvivenza più potente",
            "BAX":   "BCL-2-Associated X; pro-apoptotico; forma pori sull'OMM causando rilascio del citocromo c",
            "BAK1":  "BCL-2 Antagonist/Killer 1; pro-apoptotico residente nell'OMM; amplifica il segnale di BAX",
        }
    },
}

#dataset
def load_mitocarta(path):
    wb  = xlrd.open_workbook(path)
    sh  = wb.sheet_by_name("A Human MitoCarta3.0")
    hdr = sh.row_values(0)

    sym_col  = hdr.index("Symbol")
    desc_col = hdr.index("Description")
    loc_col  = hdr.index("MitoCarta3.0_SubMitoLocalization")
    path_col = hdr.index("MitoCarta3.0_MitoPathways")
    score_col= hdr.index("MitoCarta2.0_Score")

    rows = []
    for r in range(1, sh.nrows):
        rows.append({
            "gene":         str(sh.cell_value(r, sym_col)).strip(),
            "description":  str(sh.cell_value(r, desc_col)).strip(),
            "sublocation":  str(sh.cell_value(r, loc_col)).strip(),
            "pathways":     str(sh.cell_value(r, path_col)).strip(),
            "mc_score":     sh.cell_value(r, score_col),
        })

    sh_b = wb.sheet_by_name("B Human All Genes")
    hdr_b = sh_b.row_values(0)
    sym_b  = hdr_b.index("Symbol")
    sc_b   = hdr_b.index("MitoCarta2.0_Score")
    all_scores = {
        str(sh_b.cell_value(r, sym_b)).strip(): sh_b.cell_value(r, sc_b)
        for r in range(1, sh_b.nrows)
    }

    df_mc = pd.DataFrame(rows)
    df_mc["mc_score"] = pd.to_numeric(df_mc["mc_score"], errors="coerce").fillna(0.0)
    df_mc["gene_upper"] = df_mc["gene"].str.upper()
    confirmed = set(df_mc["gene_upper"])
    return df_mc, confirmed, all_scores


def load_surfaceome(path):
    wb  = xlrd.open_workbook(path)
    sh  = wb.sheet_by_name("11.7_Surfaceome")
    hdr = sh.row_values(1)   # riga 0 = titolo, riga 1 = intestazione

    gene_col  = hdr.index("UniProt gene")
    desc_col  = hdr.index("UniProt description")
    label_col = hdr.index("Surfaceome Label")
    src_col   = hdr.index("Surfaceome Label Source")

    rows = []
    for r in range(2, sh.nrows):
        gene = str(sh.cell_value(r, gene_col)).strip()
        if gene and gene.upper() != "NAN":
            rows.append({
                "gene":        gene,
                "description": str(sh.cell_value(r, desc_col)).strip(),
                "label":       str(sh.cell_value(r, label_col)).strip(),
                "source":      str(sh.cell_value(r, src_col)).strip(),
            })

    df_surf = pd.DataFrame(rows)
    df_surf["gene_upper"] = df_surf["gene"].str.upper()
    # deduplica (un gene può avere più isoforme)
    df_surf = df_surf.sort_values("gene_upper").drop_duplicates("gene_upper", keep="first")
    surf_set = set(df_surf["gene_upper"])
    surf_map = df_surf.set_index("gene_upper")["description"].to_dict()
    return df_surf, surf_set, surf_map

#main

def main():
    print("=" * 72)
    print("GENI ASSOCIATI ALLA DONAZIONE MITOCONDRIALE SECONDARIA")
    print("Fonti: MitoCarta 3.0  ×  PNAS Surfaceome Dataset S11.7")
    print("=" * 72)

    print("\n[1] Caricamento MitoCarta 3.0 …")
    df_mc, mc_confirmed, mc_all_scores = load_mitocarta(
        os.path.join(DATA, "Human.MitoCarta3.0.xls")
    )
    print(f"    {len(mc_confirmed)} geni mitocondriali confermati")

    print("[2] Caricamento PNAS Surfaceome …")
    df_surf, surf_set, surf_desc_map = load_surfaceome(
        os.path.join(DATA, "pnas.1808790115.sd01.xls")
    )
    print(f"    {len(surf_set)} proteine di superficie annotate")

    mc_desc_map = df_mc.set_index("gene_upper")["description"].to_dict()
    mc_loc_map  = df_mc.set_index("gene_upper")["sublocation"].to_dict()
    mc_path_map = df_mc.set_index("gene_upper")["pathways"].to_dict()

    records = []
    for module_name, info in MODULES.items():
        for gene, our_description in info["genes"].items():
            gu = gene.upper()
            in_mc   = gu in mc_confirmed
            in_surf = gu in surf_set
            records.append({
                "gene":             gene,
                "module":           module_name,
                "weight":           info["weight"],
                "category":         info["category"],
                "funzione_nostra":  our_description,
                "in_MitoCarta":     in_mc,
                "mc_score":         round(mc_all_scores.get(gu, 0.0), 2),
                "mc_sublocation":   mc_loc_map.get(gu, "—"),
                "mc_pathways":      mc_path_map.get(gu, "—"),
                "mc_description":   mc_desc_map.get(gu, "—"),
                "in_Surfaceome":    in_surf,
                "surf_description": surf_desc_map.get(gu, "—"),
            })

    df_genes = pd.DataFrame(records)

    csv_path = os.path.join(OUT, "genes_mito_donation.csv")
    df_genes.to_csv(csv_path, index=False, encoding="utf-8")
    print(f"\n[3] Tabella salvata: {csv_path}")
    print(f"    {len(df_genes)} geni totali | "
          f"{df_genes['in_MitoCarta'].sum()} in MitoCarta | "
          f"{df_genes['in_Surfaceome'].sum()} in Surfaceome")

    # report
    lines = []
    lines.append("=" * 72)
    lines.append("GENI ASSOCIATI ALLA DONAZIONE MITOCONDRIALE SECONDARIA")
    lines.append("MitoCarta 3.0  ×  PNAS Surfaceome  |  Pipeline Über Score")
    lines.append("=" * 72)
    lines.append("")
    lines.append("LEGENDA COLONNE:")
    lines.append("  [MC]   = confermato in MitoCarta 3.0 (proteina mitocondriale)")
    lines.append("  [SURF] = confermato nel PNAS Surfaceome (proteina di superficie)")
    lines.append("  Score  = punteggio probabilistico MitoCarta 2.0 (0–50+)")
    lines.append("")

    for module_name, info in MODULES.items():
        weight_str = f"w={info['weight']}" if info["weight"] else "separato"
        lines.append("─" * 72)
        lines.append(f"MODULO: {module_name}  [{weight_str}]  [{info['category'].upper()}]")
        lines.append("─" * 72)

        for gene, our_desc in info["genes"].items():
            gu = gene.upper()
            in_mc   = gu in mc_confirmed
            in_surf = gu in surf_set
            score   = mc_all_scores.get(gu, 0.0)
            loc     = mc_loc_map.get(gu, "")
            pathway = mc_path_map.get(gu, "")
            mc_d    = mc_desc_map.get(gu, "")
            surf_d  = surf_desc_map.get(gu, "")

            tags = []
            if in_mc:   tags.append(f"[MC score={score:.1f}]")
            if in_surf: tags.append("[SURF]")
            tag_str = "  " + "  ".join(tags) if tags else ""

            lines.append(f"\n  {gene}{tag_str}")
            lines.append(f"    Funzione (donazione mitocondriale):")
            import textwrap
            for l in textwrap.wrap(our_desc, 66):
                lines.append(f"      {l}")
            if mc_d and mc_d != "—":
                lines.append(f"    MitoCarta descrizione: {mc_d}")
            if loc and loc != "—":
                lines.append(f"    Localizzazione sub-mito: {loc}")
            if pathway and pathway != "—":
                top_paths = " | ".join(
                    p.strip() for p in pathway.split("|")[:3]
                )
                lines.append(f"    MitoPathways: {top_paths}")
            if surf_d and surf_d != "—":
                lines.append(f"    Surfaceome descrizione: {surf_d}")

        lines.append("")

    lines.append("=" * 72)
    lines.append("GENI MITOCARTA 3.0 NEI PATHWAYS CHIAVE PER TRASFERIMENTO")
    lines.append("(tutti i geni confermati, non solo quelli del pipeline)")
    lines.append("=" * 72)

    KEY_PATHWAYS = {
        "Mitochondrial dynamics": "Fusione, fissione, morfologia della rete mitocondriale",
        "Mitophagy":              "Degradazione selettiva di mitocondri danneggiati",
        "Protein import":         "Importazione di proteine dal citoplasma nel mitocondrio",
        "OXPHOS":                 "Catena respiratoria e sintesi di ATP",
        "Metabolism > Lipid":     "Metabolismo lipidico mitocondriale (β-ossidazione, ecc.)",
        "mtDNA":                  "Replicazione e mantenimento del DNA mitocondriale",
    }

    for kw, description in KEY_PATHWAYS.items():
        mask = df_mc["pathways"].str.contains(kw, case=False, na=False)
        subset = df_mc[mask].sort_values("mc_score", ascending=False)
        lines.append(f"\n  {kw}  —  {description}")
        lines.append(f"  ({len(subset)} geni in MitoCarta)")
        top = subset.head(10)
        for _, row in top.iterrows():
            lines.append(
                f"    {row['gene']:<12} score={row['mc_score']:5.1f}  "
                f"loc={row['sublocation']:<6}  {row['description'][:55]}"
            )

    lines.append("")
    lines.append("=" * 72)
    lines.append("PROTEINE DI SUPERFICIE (PNAS SURFACEOME) RILEVANTI AL TRASFERIMENTO")
    lines.append("=" * 72)

    SURF_KEYWORDS = [
        ("EXOSOME", "Biogenesi e rilascio di esosomi"),
        ("TETRASPANIN", "Tetraspanici organizzatori di esosomi"),
        ("INTEGRIN",    "Integrine per adesione e contatto cellula-cellula"),
        ("CHEMOKINE",   "Recettori chemochinergici per homing e migrazione"),
        ("LAMP",        "Proteine lisosomali esportate durante degranulazione"),
        ("CD63|CD9|CD81|CD151|CD82", "Tetraspanici canonici degli esosomi"),
    ]

    for kw, label in SURF_KEYWORDS:
        mask = (
            df_surf["description"].str.contains(kw, case=False, na=False) |
            df_surf["gene"].str.contains(kw, case=False, na=False)
        )
        subset = df_surf[mask]
        if subset.empty:
            continue
        lines.append(f"\n  {label}  ({len(subset)} proteine)")
        for _, row in subset.head(12).iterrows():
            lines.append(
                f"    {row['gene']:<12}  {row['description'][:58]}"
            )

    txt_path = os.path.join(OUT, "genes_mito_donation.txt")
    with open(txt_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(f"[4] Report narrativo salvato: {txt_path}")

    print()
    print("SOMMARIO PER MODULO:")
    print(f"  {'Modulo':<42} {'Geni':>5}  {'in MC':>5}  {'in Surf':>7}")
    print("  " + "-" * 64)
    for module_name, info in MODULES.items():
        genes = list(info["genes"].keys())
        n_mc   = sum(1 for g in genes if g.upper() in mc_confirmed)
        n_surf = sum(1 for g in genes if g.upper() in surf_set)
        print(f"  {module_name:<42} {len(genes):>5}  {n_mc:>5}  {n_surf:>7}")

    print()
    print("TOP GENI IN ENTRAMBI I DATABASE (MitoCarta + Surfaceome):")
    both = df_genes[df_genes["in_MitoCarta"] & df_genes["in_Surfaceome"]]
    for _, row in both.iterrows():
        print(f"  {row['gene']:<12}  [{row['module']}]")
        print(f"             MC: {row['mc_description'][:55]}")
        print(f"           Surf: {row['surf_description'][:55]}")
        print()


if __name__ == "__main__":
    main()
