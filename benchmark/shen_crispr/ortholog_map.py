"""
P4 Validation — Human-Mouse Ortholog Mapping
=============================================
Maps human gene symbols (HGNC) from IDS-SALS frozen results to mouse
gene symbols (MGI) for cross-species comparison with perturbation data.

Strategy:
  1. Primary: Biomart REST API query (Ensembl)
  2. Fallback: Simple capitalization heuristic (Human: GENE → Mouse: Gene)
  3. Cache: JSON file to avoid repeated API calls

The IDS-SALS frozen results use ENSEMBL versioned IDs (e.g., ENSG00000114541.15)
with gene symbols in the stable_high_genes.csv. The perturbation data uses
mouse gene symbols directly (e.g., Trem2, C9orf72, Lrrk2).
"""
import json
import logging
import re
import urllib.request
import urllib.error
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ── Known manual mappings (genes where simple capitalization doesn't work) ──
# These are the 18 perturbation target genes from Shen et al.
# Human HGNC symbol → Mouse MGI symbol
MANUAL_ORTHOLOGS = {
    # ALS genes
    "C9orf72": "C9orf72",      # Same symbol in mouse
    "TBK1": "Tbk1",
    "CFAP410": "Cfap410",      # a.k.a. C21orf2
    "DPP6": "Dpp6",
    # PD genes
    "LRRK2": "Lrrk2",
    "STK39": "Stk39",
    "SH3GL2": "Sh3gl2",
    # AD genes
    "TREM2": "Trem2",
    "CLU": "Clu",
    "RBFOX1": "Rbfox1",
    "NDUFAF2": "Ndufaf2",
    # Control genes
    "FASN": "Fasn",
    "FLCN": "Flcn",
    "GFAP": "Gfap",
    "OLIG2": "Olig2",
    "RRAGA": "Rraga",
    "SRF": "Srf",
    # Ribosomal / housekeeping (commonly in top φ)
    "RPL13A": "Rpl13a",
    "RPS15": "Rps15",
    "RPL13": "Rpl13",
    "RPS14": "Rps14",
    # Other common IDS-SALS genes
    "GSTK1": "Gstk1",
    "FRMD4B": "Frmd4b",
    "CSRNP3": "Csrnp3",
    "EML1": "Eml1",
    "NACAD": "Nacad",
    "CPB2": "Cpb2",
    "ENO4": "Eno4",
}


def capitalize_heuristic(human_symbol: str) -> str:
    """Simple heuristic: Human UPPER → Mouse Title case.

    Works for most protein-coding genes (e.g., TREM2 → Trem2).
    Handles special patterns:
      - Numeric suffixes (RPL13A → Rpl13a)
      - orf notation (C9orf72 → C9orf72)
      - lncRNA with -AS suffix preserved
    """
    s = human_symbol.strip()
    if not s:
        return s

    # Special patterns: C[0-9]orf[0-9]+
    m = re.match(r'^(C\d+orf)(\d+)(.*)$', s, re.IGNORECASE)
    if m:
        return m.group(1).capitalize() + m.group(2) + m.group(3).lower()

    # lncRNA with -AS suffix
    if "-AS" in s.upper():
        base, _, suffix = s.partition("-")
        return base[0].upper() + base[1:].lower() + "-" + suffix

    # Standard: first letter upper, rest lower
    return s[0].upper() + s[1:].lower()


def fetch_biomart_orthologs(
    human_symbols: List[str],
    cache_path: Optional[Path] = None,
) -> Dict[str, str]:
    """Fetch human-mouse orthologs via Ensembl BioMart REST API.

    Parameters
    ----------
    human_symbols : list of str
        Human gene symbols to query.
    cache_path : Path, optional
        If provided, cache results to this JSON file.

    Returns
    -------
    dict : {human_symbol: mouse_symbol}
    """
    if cache_path and cache_path.exists():
        logger.info(f"Loading cached orthologs from {cache_path}")
        with open(cache_path) as f:
            return json.load(f)

    logger.info(f"Querying BioMart for {len(human_symbols)} human genes...")

    # BioMart XML query
    xml_template = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1"
       uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="hgnc_symbol" value="{symbols}"/>
    <Attribute name="hgnc_symbol"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="mmusculus_homolog_associated_gene_name"/>
    <Attribute name="mmusculus_homolog_orthology_type"/>
    <Attribute name="mmusculus_homolog_orthology_confidence"/>
  </Dataset>
</Query>'''

    # BioMart has a limit on query length; batch if needed
    batch_size = 200
    result = {}

    for i in range(0, len(human_symbols), batch_size):
        batch = human_symbols[i:i + batch_size]
        symbols_str = ",".join(batch)
        xml = xml_template.format(symbols=symbols_str)

        url = "https://www.ensembl.org/biomart/martservice"
        try:
            data = urllib.parse.urlencode({"query": xml}).encode()
            req = urllib.request.Request(url, data=data)
            with urllib.request.urlopen(req, timeout=60) as resp:
                text = resp.read().decode("utf-8")

            lines = text.strip().split("\n")
            if len(lines) < 2:
                logger.warning(f"BioMart batch {i}: empty response")
                continue

            header = lines[0].split("\t")
            for line in lines[1:]:
                cols = line.split("\t")
                if len(cols) < 5:
                    continue
                human_sym = cols[0].strip()
                mouse_sym = cols[2].strip()
                orth_type = cols[3].strip()
                confidence = cols[4].strip()

                if not human_sym or not mouse_sym:
                    continue

                # Prefer one2one orthologs with high confidence
                if human_sym in result:
                    # Keep existing if it was one2one
                    if orth_type != "ortholog_one2one":
                        continue
                result[human_sym] = mouse_sym

            logger.info(f"BioMart batch {i}: {len(batch)} queried, {len(result)} mapped so far")

        except (urllib.error.URLError, TimeoutError) as e:
            logger.warning(f"BioMart query failed: {e}")
            break

    if cache_path and result:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(result, f, indent=2)
        logger.info(f"Cached {len(result)} orthologs to {cache_path}")

    return result


def build_ortholog_map(
    human_symbols: List[str],
    use_biomart: bool = True,
    cache_path: Optional[Path] = None,
) -> Dict[str, str]:
    """Build comprehensive human→mouse ortholog mapping.

    Priority:
      1. Manual curated mappings (MANUAL_ORTHOLOGS)
      2. BioMart API (if use_biomart=True)
      3. Capitalization heuristic

    Parameters
    ----------
    human_symbols : list of str
        Human gene symbols to map.
    use_biomart : bool
        Whether to query BioMart API.
    cache_path : Path, optional
        Cache file for BioMart results.

    Returns
    -------
    dict : {human_symbol: mouse_symbol}
    """
    mapping = {}
    unmapped = []

    # Step 1: Manual mappings
    for sym in human_symbols:
        if sym in MANUAL_ORTHOLOGS:
            mapping[sym] = MANUAL_ORTHOLOGS[sym]
        else:
            unmapped.append(sym)

    logger.info(f"Manual mappings: {len(mapping)}/{len(human_symbols)}")

    # Step 2: BioMart
    if use_biomart and unmapped:
        biomart = fetch_biomart_orthologs(unmapped, cache_path=cache_path)
        for sym in list(unmapped):
            if sym in biomart:
                mapping[sym] = biomart[sym]
                unmapped.remove(sym)
        logger.info(f"After BioMart: {len(mapping)}/{len(human_symbols)} "
                     f"({len(unmapped)} remaining)")

    # Step 3: Capitalization heuristic for remainder
    for sym in unmapped:
        mapping[sym] = capitalize_heuristic(sym)

    logger.info(f"Final mapping: {len(mapping)}/{len(human_symbols)} "
                 f"(heuristic fallback: {len(unmapped)})")

    return mapping


def build_reverse_map(
    ortholog_map: Dict[str, str],
) -> Dict[str, List[str]]:
    """Build mouse→human reverse mapping (one mouse gene may map to multiple human).

    Returns
    -------
    dict : {mouse_symbol: [human_symbol, ...]}
    """
    rev = {}
    for h, m in ortholog_map.items():
        rev.setdefault(m, []).append(h)
    return rev


def map_ensembl_to_symbol(
    stable_high_path: str,
) -> Dict[str, str]:
    """Extract ENSEMBL versioned ID → HGNC symbol mapping from stable_high_genes.csv.

    Returns
    -------
    dict : {ensg_versioned: symbol}
    """
    import csv
    mapping = {}
    with open(stable_high_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ensg = row.get("ensg_versioned", "").strip()
            symbol = row.get("symbol", "").strip()
            if ensg and symbol:
                mapping[ensg] = symbol
    return mapping


def fetch_mygene_ensg_to_symbol(
    ensembl_ids: List[str],
    cache_path: Optional[Path] = None,
) -> Dict[str, str]:
    """Resolve ENSEMBL gene IDs to HGNC symbols via MyGene.info REST API.

    MyGene.info is more reliable than BioMart for batch queries.

    Parameters
    ----------
    ensembl_ids : list of str
        ENSEMBL gene IDs (versioned OK — version suffix stripped before query).
    cache_path : Path, optional
        Cache file for results.

    Returns
    -------
    dict : {ensg_bare: symbol}
    """
    if cache_path and cache_path.exists():
        logger.info(f"Loading cached ENSG→symbol from {cache_path}")
        with open(cache_path) as f:
            return json.load(f)

    # Strip version suffixes
    bare_ids = list(set(eid.split(".")[0] for eid in ensembl_ids))
    logger.info(f"Querying MyGene.info for {len(bare_ids)} ENSEMBL IDs...")

    import urllib.parse

    url = "https://mygene.info/v3/query"
    batch_size = 500
    result = {}

    for i in range(0, len(bare_ids), batch_size):
        batch = bare_ids[i:i + batch_size]
        post_data = urllib.parse.urlencode({
            "q": ",".join(batch),
            "scopes": "ensembl.gene",
            "fields": "symbol",
            "species": "human",
            "size": str(len(batch)),
        }).encode()

        try:
            req = urllib.request.Request(url, data=post_data, method="POST")
            req.add_header("Content-Type", "application/x-www-form-urlencoded")
            with urllib.request.urlopen(req, timeout=60) as resp:
                items = json.loads(resp.read().decode("utf-8"))

            for item in items:
                q = item.get("query", "")
                sym = item.get("symbol", "")
                if q and sym and not item.get("notfound", False):
                    result[q] = sym

            logger.info(f"MyGene batch {i}: {len(batch)} queried, "
                         f"{len(result)} resolved so far")

        except (urllib.error.URLError, TimeoutError) as e:
            logger.warning(f"MyGene query failed: {e}")
            break

    if cache_path and result:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(result, f, indent=2)
        logger.info(f"Cached {len(result)} ENSG→symbol to {cache_path}")

    return result


def resolve_ensg_to_symbol(
    ensembl_ids: List[str],
    frozen_results_dir: str,
    gene_annotation_path: Optional[str] = None,
    biomart_cache_path: Optional[Path] = None,
) -> Dict[str, str]:
    """Resolve ENSEMBL versioned IDs to HGNC symbols.

    Priority:
      1. stable_high_genes.csv (has 135 genes with symbols)
      2. gene_annotation_utf8.csv (if available on analysis machine)
      3. BioMart API (for remaining)

    Parameters
    ----------
    ensembl_ids : list of str
        ENSEMBL versioned IDs (e.g., "ENSG00000114541.15").
    frozen_results_dir : str
        Path to frozen results directory.
    gene_annotation_path : str, optional
        Path to gene_annotation_utf8.csv (from analysis machine resources/).
    biomart_cache_path : Path, optional
        Cache file for BioMart results.

    Returns
    -------
    dict : {ensg_versioned: symbol}
    """
    import csv

    mapping = {}

    # Step 1: stable_high_genes.csv
    stable_path = (Path(frozen_results_dir) / "track_b"
                   / "sals_upstream_gene_list" / "stable_high_genes.csv")
    if stable_path.exists():
        with open(stable_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                ensg = row.get("ensg_versioned", "").strip()
                symbol = row.get("symbol", "").strip()
                if ensg and symbol:
                    mapping[ensg] = symbol
    logger.info(f"Step 1 (stable_high): {len(mapping)} symbols")

    # Step 2: gene_annotation_utf8.csv (if available)
    if gene_annotation_path and Path(gene_annotation_path).exists():
        import pandas as pd
        try:
            annot = pd.read_csv(gene_annotation_path)
            # Expected columns: gene_id (bare ENSG), gene_name (symbol)
            if "gene_id" in annot.columns and "gene_name" in annot.columns:
                bare_to_sym = dict(zip(annot["gene_id"], annot["gene_name"]))
                for eid in ensembl_ids:
                    if eid not in mapping:
                        bare = eid.split(".")[0]
                        if bare in bare_to_sym:
                            mapping[eid] = bare_to_sym[bare]
        except Exception as e:
            logger.warning(f"gene_annotation load failed: {e}")
        logger.info(f"Step 2 (annotation): {len(mapping)} symbols")

    # Step 3: MyGene.info API for remainder (more reliable than BioMart)
    unmapped_ids = [eid for eid in ensembl_ids if eid not in mapping]
    if unmapped_ids:
        mygene = fetch_mygene_ensg_to_symbol(unmapped_ids, cache_path=biomart_cache_path)
        for eid in unmapped_ids:
            bare = eid.split(".")[0]
            if bare in mygene:
                mapping[eid] = mygene[bare]
        logger.info(f"Step 3 (MyGene.info): {len(mapping)} symbols total")

    still_unmapped = [eid for eid in ensembl_ids if eid not in mapping]
    if still_unmapped:
        logger.warning(f"{len(still_unmapped)} ENSEMBL IDs could not be resolved to symbols")

    return mapping


# ── Convenience: Build everything from frozen results ──

def load_frozen_gene_rankings(
    frozen_results_dir: str,
    source: str = "stable_high",
) -> Tuple[Dict[str, float], Dict[str, str], List[str]]:
    """Load gene-level φ rankings from frozen IDS-SALS results.

    Parameters
    ----------
    frozen_results_dir : str
        Path to results/ directory in frozen analysis.
    source : str
        "stable_high" - only 135 bootstrap-stable high-φ genes (default)
        "rank_shift"  - all 922 genes from rank_shift_per_gene.csv

    Returns
    -------
    phi_by_symbol : dict
        {HGNC_symbol: phi_score}
    ensg_to_symbol : dict
        {ENSG_versioned: HGNC_symbol}
    ranked_symbols : list
        Symbols sorted by φ descending
    """
    import csv

    if source == "stable_high":
        stable_path = (Path(frozen_results_dir) / "track_b"
                       / "sals_upstream_gene_list" / "stable_high_genes.csv")
        if not stable_path.exists():
            raise FileNotFoundError(f"stable_high_genes.csv not found at {stable_path}")

        phi_by_symbol = {}
        ensg_to_symbol = {}
        ranked = []

        with open(stable_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                ensg = row.get("ensg_versioned", "").strip()
                symbol = row.get("symbol", "").strip()
                phi = float(row.get("phi", 0))
                if ensg and symbol:
                    ensg_to_symbol[ensg] = symbol
                    phi_by_symbol[symbol] = phi
                    ranked.append((symbol, phi))

        ranked.sort(key=lambda x: -x[1])
        ranked_symbols = [s for s, _ in ranked]
        logger.info(f"Loaded {len(phi_by_symbol)} gene φ scores (stable_high)")
        return phi_by_symbol, ensg_to_symbol, ranked_symbols

    elif source == "rank_shift":
        # Load full 922-gene φ from rank_shift_per_gene.csv
        # Uses phi_disease (SALS condition) as the primary φ
        rs_path = (Path(frozen_results_dir) / "track_b" / "laneB"
                   / "rank_shift_oligo_edge_weight" / "rank_shift_per_gene.csv")
        if not rs_path.exists():
            raise FileNotFoundError(f"rank_shift_per_gene.csv not found at {rs_path}")

        phi_by_ensg = {}
        ensg_list = []
        with open(rs_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                ensg = row["gene"].strip()
                phi_disease = float(row["phi_disease"])
                phi_by_ensg[ensg] = phi_disease
                ensg_list.append(ensg)

        logger.info(f"Loaded {len(phi_by_ensg)} genes from rank_shift (φ_disease)")
        # Note: ENSG→symbol resolution must be done separately
        # Return with ENSG as keys (caller must resolve symbols)
        return phi_by_ensg, {}, ensg_list  # ensg_to_symbol empty, ranked is ensg list

    else:
        raise ValueError(f"Unknown source: {source}")


def load_rank_shift_per_gene(
    frozen_results_dir: str,
    celltype: str = "oligo",
    flow_mode: str = "edge_weight",
) -> Dict[str, Dict]:
    """Load per-gene rank shift results for a specific celltype.

    Returns
    -------
    dict : {ENSG_versioned: {phi_disease, phi_control, rank_shift, ...}}
    """
    import csv

    rs_dir = (Path(frozen_results_dir) / "track_b" / "laneB"
              / f"rank_shift_{celltype}_{flow_mode}")
    rs_path = rs_dir / "rank_shift_per_gene.csv"

    if not rs_path.exists():
        raise FileNotFoundError(f"rank_shift_per_gene.csv not found at {rs_path}")

    result = {}
    with open(rs_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene = row["gene"].strip()
            result[gene] = {
                "phi_disease": float(row["phi_disease"]),
                "phi_control": float(row["phi_control"]),
                "rank_disease": float(row["rank_disease"]),
                "rank_control": float(row["rank_control"]),
                "rank_shift": float(row["rank_shift"]),
                "p_value": float(row["p_value"]),
                "classification": row["classification"].strip(),
            }

    logger.info(f"Loaded {len(result)} rank shift entries from {rs_path}")
    return result
