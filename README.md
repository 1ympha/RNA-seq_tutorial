# RNA-seq_tutorial
# RNA-seq Analysis Pipeline (fastp → SortMeRNA → HISAT2 → StringTie)

This repository provides a **one-stop, end-to-end RNA-seq analysis workflow** for paired-end Illumina data, from raw FASTQ files to gene-level TPM and count matrices.

- Species: *Mus musculus*
- Genome version: **GRCm39**
- Sequencing type: Paired-end RNA-seq
- Audience: Beginners to intermediate users

---

## 0. Environment Setup

Activate conda environment:

```bash
source /public/home/shiyunfeng/software/anaconda3/bin/activate
```

Required software:

- fastp
- fastqc
- multiqc
- sortmerna
- hisat2
- samtools
- stringtie
- python ≥ 3.7 (pandas, tqdm)

---

## 1. Raw Data Quality Control (fastp)

### Purpose
- Adapter trimming
- Low-quality read filtering
- Maintain paired-end consistency

```bash
#!/bin/bash
set -euo pipefail

INDIR="path/to/rawdata"
OUTDIR="path/to/cleanData"

mkdir -p "$OUTDIR"

for r1 in "$INDIR"/*.R1.fq.gz; do
    base=$(basename "$r1" .R1.fq.gz)
    r2="$INDIR/${base}.R2.fq.gz"

    fastp \
        -i "$r1" -I "$r2" \
        -o "$OUTDIR/${base}.R1.clean.fq.gz" \
        -O "$OUTDIR/${base}.R2.clean.fq.gz" \
        --thread 10
done
```

---

## 2. Quality Assessment (FastQC + MultiQC)

```bash
#!/bin/bash
set -euo pipefail

CLEANDIR="path/to/cleanData"
OUTDIR="path/to/qc_multiqc_report"

mkdir -p "$OUTDIR"

for r1 in "$CLEANDIR"/*.R1.clean.fq.gz; do
    base=$(basename "$r1" .R1.clean.fq.gz)
    r2="$CLEANDIR/${base}.R2.clean.fq.gz"

    fastqc "$r1" "$r2" -t 10 -o "$OUTDIR"
done

multiqc "$OUTDIR" -o "$OUTDIR" -n multiqc
```

---

## 3. rRNA Removal (SortMeRNA)

```bash
ls cleanData/*.R1.clean.fq.gz > file.list
split -l 40 file.list sl_
```

```bash
#!/bin/bash
set -euo pipefail

SOFTWARE_DIR="/public/software/sortmerna-4.0.0/data/rRNA_databases"
FILELIST="$1"

for r1 in $(cat "$FILELIST"); do
    r2=${r1/R1./R2.}
    sample=$(basename "$r1" .R1.clean.fq.gz)

    sortmerna \
        --threads 20 \
        --ref ${SOFTWARE_DIR}/silva-bac-16s-id90.fasta \
        --ref ${SOFTWARE_DIR}/silva-bac-23s-id98.fasta \
        --ref ${SOFTWARE_DIR}/silva-arc-16s-id95.fasta \
        --ref ${SOFTWARE_DIR}/silva-arc-23s-id98.fasta \
        --ref ${SOFTWARE_DIR}/silva-euk-18s-id95.fasta \
        --ref ${SOFTWARE_DIR}/silva-euk-28s-id98.fasta \
        --ref ${SOFTWARE_DIR}/rfam-5s-database-id98.fasta \
        --ref ${SOFTWARE_DIR}/rfam-5.8s-database-id98.fasta \
        --reads "$r1" --reads "$r2" \
        --paired_in --fastx --out2 \
        --workdir sortmerna/run/"$sample" \
        --aligned sortmerna/output_rRNA/"$sample" \
        --other sortmerna/sortmerna_clean/"$sample"
done
```

---

## 4. Alignment (HISAT2)

```bash
#!/bin/bash
set -euo pipefail

INDIR="path/to/sortmerna_clean"
OUTDIR="path/to/alignment"
IDX="/public/database/Genome_File/Mus_musculus.GRCm39.hisat2.index/genome"

mkdir -p "$OUTDIR"

for r1 in "$INDIR"/*.R1.clean.fq.gz; do
    sample=$(basename "$r1" .R1.clean.fq.gz)
    r2="$INDIR/${sample}.R2.clean.fq.gz"

    hisat2 \
        --dta -p 10 -t \
        -x "$IDX" \
        -1 "$r1" -2 "$r2" \
        -S "$OUTDIR/${sample}.sam"
done
```

---

## 5. Mapping Statistics

```bash
cd path/to/alignment

for sam in *.sam; do
    samtools flagstat "$sam" > "${sam%.sam}.stats"
done
```

```bash
echo -e "Sample\tMapping_Rate" > mapping_rates_summary.tsv

for stat in *.stats; do
    sample=${stat%.stats}
    rate=$(grep "mapped (" "$stat" | awk -F '[()%]' '{print $2}')
    echo -e "$sample\t$rate" >> mapping_rates_summary.tsv
done
```

---

## 6. SAM to BAM and Sorting

```bash
cd path/to/alignment

for sam in *.sam; do
    base=${sam%.sam}
    samtools view -@ 15 -bhS -q 30 "$sam" > "${base}.uniq.bam"
    samtools sort -@ 15 -o "${base}.uniq.sorted.bam" "${base}.uniq.bam"
done
```

---

## 7. Quantification (StringTie)

```bash
#!/bin/bash
set -euo pipefail

BAMDIR="path/to/alignment"
GTF="/public/database/Genome_File/Mouse_genomeFile/Mus_musculus.GRCm39.114.gtf"
OUTDIR="$BAMDIR/gtf_2"

mkdir -p "$OUTDIR"
cd "$BAMDIR"

for bam in *.uniq.sorted.bam; do
    sample=${bam%.uniq.sorted.bam}

    stringtie "$bam" \
        -p 15 \
        -e -G "$GTF" \
        -o "$OUTDIR/${sample}.gtf" \
        -A "$OUTDIR/${sample}.gene_abundance.txt"
done
```

---

## 8. Merge Count Matrix (prepDE.py)

```bash
cd gtf_2
ls *.gtf | awk -F'.gtf' '{print $1"\t"$0}' > sample_list.txt
```

```bash
python3 prepDE.py3 \
    -i sample_list.txt \
    -g gene_count_matrix.csv \
    -t transcript_count_matrix.csv
```

---

## 9. Merge TPM Matrix

```python
import os
import pandas as pd
from tqdm import tqdm

input_dir = "path/to/alignment/gtf_2"
dfs = []

for f in tqdm(os.listdir(input_dir)):
    if f.endswith(".gene_abundance.txt"):
        sample = f.split(".")[0]
        df = pd.read_csv(os.path.join(input_dir, f), sep="\\t")
        df = df.drop_duplicates("Gene ID").set_index("Gene ID")[["TPM"]]
        df.columns = [f"{sample}_TPM"]
        dfs.append(df)

result = pd.concat(dfs, axis=1)
result.to_csv("TPM_All_Samples.txt", sep="\\t")
```

---

## Final Outputs

- `multiqc.html`
- `mapping_rates_summary.tsv`
- `gene_count_matrix.csv`
- `TPM_All_Samples.txt`
