# 🧬 FASTA/FASTQ-Based Variant Annotator App

This Streamlit app allows you to upload a FASTA or FASTQ file, align it to the GRCh38 reference genome, call variants, annotate them using multiple genomic databases (NCBI, MyVariant.info, Ensembl, UCSC), and retrieve disease + pharmaceutical information with visualizations.

## ✅ Features
- Input: FASTA/FASTQ
- Alignment: BWA to GRCh38
- Variant Calling: samtools + bcftools
- Annotation via rsID: NCBI → MyVariant → Ensembl → UCSC fallback
- Disease mapping
- Clinical trial lookup
- Visualization and CSV download

## 🔧 Installation
```bash
pip install -r requirements.txt
conda install -c bioconda bwa samtools bcftools
```

Download GRCh38 reference genome:
```bash
mkdir data
wget -O data/hg38.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz"
gunzip data/hg38.fa.gz
```

## 🚀 Run the App
```bash
streamlit run app.py
```

## 🌐 Streamlit Cloud
This app can be deployed to [Streamlit Cloud](https://streamlit.io/cloud) by uploading:
- `app.py`
- `requirements.txt`
- `README.md`

You must skip alignment/variant-calling or use pre-made `.vcf` files, as Streamlit Cloud can't run system-level binaries (bwa, samtools).

## 📁 File Structure
```
fasta-annotator-app/
├── app.py
├── requirements.txt
├── .gitignore
├── README.md
├── data/
│   └── hg38.fa
├── tmp/
└── ... (output files)
```

## 🧪 Example Test
Upload a short `.fasta` file containing one chromosome region with known variant.
Check that annotations appear with links and clinical trial data.
