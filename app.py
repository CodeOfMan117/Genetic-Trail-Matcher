# Full FASTA/FASTQ-Based Variant Annotator App (with local variant calling)

import streamlit as st
import os
import subprocess
import pandas as pd
import plotly.express as px
import requests

# Configure app
st.set_page_config(page_title="FASTA Variant Annotator", layout="wide")
st.title("🧬 Genomic Variant Annotator from FASTA/FASTQ")

# File upload
uploaded_file = st.file_uploader("Upload FASTA or FASTQ file", type=["fasta", "fa", "fastq"])

# Output directory
os.makedirs("tmp", exist_ok=True)

# Functions for pipeline steps

def save_uploaded_file(uploaded_file):
    filepath = os.path.join("tmp", uploaded_file.name)
    with open(filepath, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return filepath

def run_alignment(input_path):
    output_sam = input_path + ".sam"
    ref = "data/hg38.fa"
    cmd = ["bwa", "mem", ref, input_path]
    with open(output_sam, "w") as out:
        subprocess.run(cmd, stdout=out)
    return output_sam

def run_variant_calling(sam_file):
    bam_file = sam_file.replace(".sam", ".bam")
    sorted_bam = bam_file.replace(".bam", ".sorted.bam")
    vcf_file = bam_file.replace(".bam", ".vcf")

    subprocess.run(["samtools", "view", "-bS", sam_file, "-o", bam_file])
    subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam])
    subprocess.run(["samtools", "index", sorted_bam])
    subprocess.run(["bcftools", "mpileup", "-f", "data/hg38.fa", sorted_bam, "-Ou"], stdout=open("tmp/tmp.bcf", "wb"))
    subprocess.run(["bcftools", "call", "-mv", "-Ov", "-o", vcf_file], input=open("tmp/tmp.bcf", "rb").read())

    return vcf_file

def extract_variants(vcf_path):
    variants = []
    with open(vcf_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom, pos, rsid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            variants.append({"chrom": chrom, "pos": pos, "rsid": rsid, "ref": ref, "alt": alt})
    return variants

def annotate_variant(variant):
    rsid = variant["rsid"]
    chrom = variant["chrom"]
    pos = variant["pos"]
    ref = variant["ref"]
    alt = variant["alt"]

    # NCBI fallback
    try:
        ncbi_url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rsid.replace('rs', '')}"
        r = requests.get(ncbi_url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            gene = data.get("primary_snapshot_data", {}).get("placements_with_allele", [{}])[0].get("alleles", [{}])[0].get("allele", {}).get("spdi", {}).get("variant_base", 'NA')
            return {
                "chr": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene or "NA",
                "clinical_significance": "NA",
                "condition": "From NCBI",
                "link": ncbi_url,
                "source": "NCBI"
            }
    except:
        pass

    # MyVariant fallback
    try:
        url_mv = f"https://myvariant.info/v1/variant/{rsid}"
        r = requests.get(url_mv, timeout=10)
        if r.status_code == 200:
            data = r.json()
            gene = data.get('gene', {}).get('symbol', 'NA')
            clinical = data.get('clinvar', {})
            return {
                "chr": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene,
                "clinical_significance": clinical.get('clinical_significance', 'Not available'),
                "condition": clinical.get('trait', ['Unknown'])[0] if isinstance(clinical.get('trait', []), list) else 'Unknown',
                "link": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinical.get('rcv', [{}])[0].get('accession', '')}" if 'rcv' in clinical else '',
                "source": "MyVariant.info"
            }
    except:
        pass

    # Ensembl fallback
    try:
        ens_url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}/{ref}/{alt}?content-type=application/json"
        r = requests.get(ens_url, timeout=10)
        if r.status_code == 200:
            ens_data = r.json()
            gene = ens_data[0]['transcript_consequences'][0].get('gene_symbol', 'NA') if ens_data else 'NA'
            return {
                "chr": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene,
                "clinical_significance": "NA",
                "condition": "From Ensembl",
                "link": ens_url,
                "source": "Ensembl"
            }
    except:
        pass

    # UCSC fallback
    return {
        "chr": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "gene": "Unknown",
        "clinical_significance": "Not available",
        "condition": "No data",
        "link": f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr{chrom}:{pos}",
        "source": "UCSC"
    }

# Start pipeline
if uploaded_file is not None:
    path = save_uploaded_file(uploaded_file)
    st.write("✅ File saved. Running alignment...")
    sam = run_alignment(path)
    vcf = run_variant_calling(sam)
    variants = extract_variants(vcf)

    st.success(f"Found {len(variants)} variants. Annotating...")
    annots = [annotate_variant(v) for v in variants]
    df = pd.DataFrame(annots)

    st.subheader("📋 Annotated Variants")
    st.dataframe(df)

    fig = px.scatter(
        df,
        x="pos",
        y=[0]*len(df),
        hover_data=["gene", "condition", "clinical_significance"],
        title="🧬 Variant Positions on Chromosome",
        labels={"pos": "Genomic Position"}
    )
    fig.update_traces(marker=dict(size=10, color="#636EFA"))
    fig.update_layout(yaxis=dict(showticklabels=False), showlegend=False)
    st.plotly_chart(fig, use_container_width=True)

    st.download_button(
        label="📥 Download Annotated CSV",
        data=df.to_csv(index=False).encode("utf-8"),
        file_name="annotated_fasta_variants.csv",
        mime="text/csv"
    )
