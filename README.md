# Pipeline de Análisis Transcriptómico RNA-Seq con Nextflow DSL2 y Despliegue en AWS Cloud

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-23aa62.svg?logo=nextflow&logoColor=white)](https://www.nextflow.io/)
[![Language-Bash](https://img.shields.io/badge/language-Bash-4EAA25.svg?logo=gnu-bash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Language-R](https://img.shields.io/badge/language-R_4.4-276DC3.svg?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status: Complete](https://img.shields.io/badge/Status-100%25_Executed-brightgreen.svg)]()

> **Bioinformática Traslacional & Cloud Computing:**  
> Implementación reproducible, modular y automatizada de un pipeline de RNA-Seq de extremo a extremo: adquisición de datos crudos (SRA), control de calidad (FastQC), filtrado y recorte de adaptadores (Trimmomatic), alineamiento guiado por empalmes (HISAT2), cuantificación génica (featureCounts) y análisis de expresión génica diferencial (DESeq2). El workflow está desarrollado bajo **Nextflow DSL2** y desacoplado arquitecturalmente para su orquestación híbrida (local vía Conda/Docker y elástica en **AWS Batch + Amazon S3**).

---

## 📑 Tabla de Contenidos
1. [Resumen del Proyecto y Fundamento Científico](#-resumen-del-proyecto-y-fundamento-cientfico)
2. [Dataset Biológico y Diseño Experimental](#-dataset-biolgico-y-diseo-experimental)
3. [Estrategia de Optimización de Almacenamiento e I/O](#-estrategia-de-optimizacin-de-almacenamiento-e-io)
4. [Arquitectura del Pipeline y Grafo de Ejecución (DAG)](#-arquitectura-del-pipeline-y-grafo-de-ejecucin-dag)
5. [Estructura del Repositorio y Descripción de Componentes](#-estructura-del-repositorio-y-descripcin-de-componentes)
6. [Resultados Bioestadísticos y Métricas Analíticas](#-resultados-bioestadsticos-y-mtricas-analticas)
7. [Guía de Ejecución Local (Conda / Docker)](#-gua-de-ejecucin-local-conda--docker)
8. [Orquestación y Despliegue en AWS Cloud (AWS Batch + S3)](#-orquestacin-y-despliegue-en-aws-cloud-aws-batch--s3)
9. [Aplicaciones en la Práctica Biomédica](#-aplicaciones-en-la-prctica-biomdica)

---

## 🧬 Resumen del Proyecto y Fundamento Científico

El análisis transcriptómico mediante **secuenciación masiva de ARN (RNA-Seq)** permite cuantificar perfiles globales de expresión génica e identificar alteraciones mecanísticas en condiciones patológicas o tras perturbaciones farmacológicas y genéticas.

Este proyecto resuelve los dos desafíos fundamentales en el procesamiento bioinformático de RNA-Seq:
1. **Reproducibilidad y escalabilidad:** Construcción de un flujo de trabajo modular con **Nextflow DSL2**, donde cada paso analítico opera en entornos aislados con control estricto de dependencias.
2. **Eficiencia en infraestructura:** Implementación de perfiles de ejecución intercambiables que permiten validar el pipeline en entornos locales limitados en hardware y escalarlo de forma nativa a la nube (**AWS Batch**) sin refactorizar una sola línea de código del pipeline.

---

## 🐭 Dataset Biológico y Diseño Experimental

El dataset analizado corresponde al estudio público del NCBI SRA **BioProject [PRJNA1533324](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1533324)**, centrado en modelos neurobiológicos de exposición prenatal al alcohol (*prenatal alcohol exposure*):

- **Organismo:** *Mus musculus* (cepa endogámica C57BL/6J).
- **Tejido:** Corteza cerebral embrionaria en estadio E18 (*embryonic day 18*).
- **Plataforma de Secuenciación:** Illumina NovaSeq 6000 (*Paired-End*, 2 × 150 pb, química con orientación *reverse-stranded*).
- **Muestras y Diseño Experimental ($n=6$):**
  - **Control Negativo ($n=3$):** `SRR40794928`, `SRR40794929`, `SRR40794930`.
  - **Inhibición de miR-150-5p ($n=3$):** `SRR40794931`, `SRR40794932`, `SRR40794933`.

**Hipótesis experimental:**  
Evaluar los cambios transcriptómicos en la corteza cerebral derivados de la perturbación del microARN regulador **miR-150-5p** frente al grupo control.

---

## 💾 Estrategia de Optimización de Almacenamiento e I/O

El conjunto crudo de lecturas de secuenciación supera los **70 GB** en formato no comprimido. En entornos computacionales restringidos, descargar simultáneamente todos los archivos colapsaría el almacenamiento antes de iniciar el procesamiento.

Para solventar esta limitación física de almacenamiento sin comprometer la integridad analítica, se desarrolló una estrategia de **adquisición y submuestreo secuencial en bucle cerrado** (`download_and_subsample.sh`):

```mermaid
flowchart LR
    A["Entrada: SRR"] --> B["fasterq-dump\n(descarga multi-hilo)"]
    B --> C["seqtk sample -s42\n(1.000.000 read pairs)"]
    C --> D["Compresión gzip\n(R1 + R2)"]
    D --> E["Purga inmediata de\narchivos crudos (>15 GB)"]
    E --> F["Siguiente muestra"]
```

- **Mapeo individual:** Cada muestra se descarga individualmente con `fasterq-dump --split-files`.
- **Submuestreo estocástico reproducible:** Mediante `seqtk sample -s42`, se extraen de manera aleatoria y reproducible exactamente **1.000.000 de pares de lecturas emparejadas** por muestra (4.000.000 de líneas FASTQ por archivo).
- **Compresión inmediata:** Generación de archivos `*.fastq.gz` (~75 MB por archivo) y eliminación completa del archivo crudo antes de iniciar la descarga de la siguiente muestra.

> **Impacto computacional:** Se redujo la huella de disco de >70 GB a **<1 GB**, permitiendo ejecutar el pipeline completo de manera local y determinista.

---

## 🏗 Arquitectura del Pipeline y Grafo de Ejecución (DAG)

El flujo de procesamiento sigue una arquitectura modular en **Nextflow DSL2**:

```
                       ┌─────────────────────────┐
                       │    samplesheet.csv      │
                       └────────────┬────────────┘
                                    │
                       [Canal: meta, [reads R1, R2]]
                                    │
                    ┌───────────────┴───────────────┐
                    ▼                               ▼
        ┌───────────────────────┐       ┌───────────────────────┐
        │      FASTQC_RAW       │       │      TRIMMOMATIC      │
        │ Control de calidad de │       │ Eliminación adapters  │
        │   lecturas crudas     │       │     TruSeq3-PE-2      │
        └───────────────────────┘       └───────────┬───────────┘
                                                    │
                                      [Lecturas procesadas R1, R2]
                                                    │
                              ┌─────────────────────┴─────────────────────┐
                              ▼                                           ▼
                  ┌───────────────────────┐                   ┌───────────────────────┐
                  │    FASTQC_TRIMMED     │                   │     HISAT2_ALIGN      │
                  │ Control de calidad de │                   │  Alineamiento a GRCm39│
                  │  lecturas recortadas  │                   │    + SAMtools sort    │
                  └───────────────────────┘                   └───────────┬───────────┘
                                                                          │
                                                                [Archivos BAM ordenados]
                                                                          │
                                                                          ▼
                                                              ┌───────────────────────┐
                                                              │     FEATURECOUNTS     │
                                                              │ Cuantificación génica │
                                                              │  (Anotación GENCODE)  │
                                                              └───────────┬───────────┘
                                                                          │
                                                                [counts.txt unificado]
                                                                          │
                                                                          ▼
                                                              ┌───────────────────────┐
                                                              │        DESEQ2         │
                                                              │ Expresión diferencial │
                                                              │  (apeglm shrinkage)   │
                                                              └───────────┬───────────┘
                                                                          │
                                           ┌──────────────────────────────┴──────────────────────────────┐
                                           ▼                                                             ▼
                               ┌───────────────────────┐                                     ┌───────────────────────┐
                               │     Volcano Plot      │                                     │       PCA Plot        │
                               │  deseq2_results.csv   │                                     │ normalized_counts.csv │
                               └───────────────────────┘                                     └───────────────────────┘
```

---

## 📂 Estructura del Repositorio y Descripción de Componentes

```
Pipeline-AWS-Nextflow/
├── main.nf                         # Orquestador del flujo Nextflow DSL2; define canales y orden de ejecución
├── nextflow.config                 # Archivo de configuración maestro (parámetros, recursos y perfiles de ejecución)
├── environment.yml                 # Especificación del entorno Conda (FastQC, Trimmomatic, HISAT2, Subread, R/Bioconductor)
├── .gitignore                      # Reglas de exclusión de Git (excluye genoma crudo de 4GB y BAMs de >100MB)
├── TruSeq3-PE-2.fa                 # Secuencias FASTA de adaptadores Illumina TruSeq para recorte por Trimmomatic
├── SRR_Acc_List.txt                # Identificadores de acceso SRA analizados
├── SraRunTable.csv                 # Metadatos experimentales descargados de NCBI SRA
│
├── modules/                        # Procesos aislados y modulares de Nextflow
│   ├── fastqc.nf                   # Módulo de control de calidad; genera métricas de distribución de calidad Phred y GC
│   ├── trimmomatic.nf              # Módulo de recorte por ventana deslizante y eliminación de adaptadores Illumina
│   ├── hisat2.nf                   # Módulo de alineamiento guiado por empalmes acoplado a SAMtools (BAM ordenado + índice BAI)
│   ├── featurecounts.nf            # Módulo de cuantificación reversa por pares de lecturas a nivel de exón/gen
│   └── deseq2.nf                   # Envoltorio del proceso de análisis estadístico en R
│
├── scripts/
│   └── deseq2_analysis.R           # Script estadístico en R: normalización, VST, PCA, apeglm shrinkage y Volcano plot
│
├── data/
│   └── samplesheet.csv             # Manifiesto estructurado (sample_id, fastq_1, fastq_2, condition)
│
└── results/                        # Entregables analíticos generados por el pipeline
    ├── deseq2/
    │   ├── deseq2_results.csv      # Tabla de resultados de expresión diferencial con contrastes y valores padj
    │   ├── normalized_counts.csv   # Matriz de conteos normalizados mediante el método median-of-ratios
    │   ├── volcano_plot.pdf        # Gráfico de volcán (Log2 Fold Change vs -Log10 p-value)
    │   └── pca_plot.pdf            # Análisis de Componentes Principales con estabilización de varianza (VST)
    ├── featurecounts/
    │   ├── counts.txt              # Matriz global de conteos crudos de todas las muestras (78.239 genes)
    │   └── counts.txt.summary      # Estadísticas de asignación de lecturas (tasa de asignación exónica vs intrónica/intergénica)
    ├── fastqc/                     # Reportes interactivos de calidad (HTML y ZIP) pre y post-procesamiento
    ├── hisat2/                     # Logs de alineamiento con porcentajes de mapeo concordante único
    ├── trimmomatic/                # Registros de filtrado con conteo de lecturas emparejadas supervivientes
    ├── pipeline_report.html        # Reporte Nextflow de consumo de CPU, memoria y tiempos de ejecución
    ├── timeline.html               # Gráfico temporal interactivo de ejecución de tareas
    └── trace.txt                   # Registro tabular detallado con métricas de cada proceso
```

---

## 📊 Resultados Bioestadísticos y Métricas Analíticas

### 1. Cuantificación y Cobertura Genómica
- **Referencia genómica:** *Mus musculus* GRCm39 primary assembly (2.6 GB) con anotación de transcritos **GENCODE vM36** (1.1 GB).
- **Alineamiento (HISAT2):** Las lecturas recortadas se mapearon de forma concordante con strandness invertido (`--rna-strandness RF`).
- **Asignación de lecturas (featureCounts):** Se evaluaron **78.239 genes**. Tras aplicar el filtro de conteo mínimo ($\ge 10$ lecturas acumuladas entre todas las réplicas) para eliminar ruido técnico estocástico, se retuvieron **13.668 genes informativos** para modelado estadístico.

### 2. Expresión Génica Diferencial (DESeq2)
El análisis implementado en [`scripts/deseq2_analysis.R`](file:///scripts/deseq2_analysis.R) ejecutó:
- **Normalización:** Estimación de factores de escala por el método de la mediana de las razones (*median-of-ratios*).
- **Modelado:** Estimación de la dispersión de conteos bajo el modelo lineal generalizado con distribución Binomial Negativa:
  $$K_{ij} \sim \text{NB}(\mu_{ij}, \alpha_i)$$
- **Contracción empírica de tasas de cambio:** Aplicación del estimador bayesiano `apeglm` sobre los $\log_2(\text{Fold Change})$ para reducir el sesgo de genes con conteos bajos.
- **Gráficos diagnósticos:**
  - **PCA:** Reducción de dimensionalidad sobre datos transformados mediante VST (*Variance Stabilizing Transformation*).
  - **Volcano plot:** Representación de significancia estadística frente a magnitud de cambio ($p_{\text{adj}} < 0.05$, $|\log_2\text{FC}| \ge 1.0$).

---

## 🚀 Guía de Ejecución Local (Conda / Docker)

### Prerrequisitos
- Sistema operativo Linux o **Windows Subsystem for Linux (WSL2 Ubuntu)**.
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) o Mamba instalado.

### 1. Clonar el repositorio
```bash
git clone https://github.com/menmati/Pipeline-AWS-Nextflow.git
cd Pipeline-AWS-Nextflow
```

### 2. Configurar el entorno bioinformático
```bash
conda env create -f environment.yml
conda activate rnaseq
```

### 3. Ejecutar el pipeline
```bash
nextflow run main.nf -profile local
```

Si deseas reanudar una ejecución interrumpida reutilizando el caché computacional:
```bash
nextflow run main.nf -profile local -resume
```

---

## ☁️ Orquestación y Despliegue en AWS Cloud (AWS Batch + S3)

> 💡 **Transparencia en infraestructura y costes:**  
> **No es necesario crear una cuenta de AWS ni incurrir en costes de facturación.**  
> Este repositorio incluye la arquitectura de nube completamente parametrizada y desacoplada mediante configuración declarativa. Permite comprender y auditar el despliegue a gran escala sin necesidad de provisionar recursos activos.

### Arquitectura de Escalado en AWS
Cuando los volúmenes de datos crecen a cientos de gigabytes o terabytes, Nextflow permite trasladar la carga computacional desde la máquina local hacia la infraestructura de Amazon Web Services sin modificar la lógica bioinformática:

```
[Estación de Trabajo / CLI / Seqera Platform]
                      │
                      │ nextflow run main.nf -profile aws
                      ▼
             ┌─────────────────┐
             │    AWS Batch    │ ◄─── Orquestador de trabajos (Job Queue)
             └────────┬────────┘
                      │
       ┌──────────────┼──────────────┐
       ▼              ▼              ▼
  ┌─────────┐    ┌─────────┐    ┌─────────┐
  │ EC2 (1) │    │ EC2 (2) │    │ EC2 (N) │ ◄── Instancias aprovisionadas elásticamente
  │ FastQC  │    │ HISAT2  │    │ DESeq2  │     con contenedores Docker de Quay/ECR
  └─────────┘    └─────────┘    └─────────┘
       │              │              │
       └──────────────┼──────────────┘
                      ▼
             ┌─────────────────┐
             │    Amazon S3    │ ◄─── Object Storage compartido
             │ s3://bucket/... │     (FASTQs, genoma, índices y resultados)
             └─────────────────┘
```

### Parametrización en [`nextflow.config`](file:///nextflow.config)

El perfil `aws` desacopla los endpoints locales y los enlaza a servicios gestionados en la nube:

```groovy
profiles {
    aws {
        process.executor  = 'awsbatch'
        process.queue     = 'rnaseq-batch-queue'
        process.container = 'quay.io/biocontainers/mulled-v2-...'

        aws {
            region = 'us-east-1'
            batch {
                cliPath = '/home/ec2-user/miniconda3/bin/aws'
                volumes = '/tmp'
            }
        }

        params {
            samplesheet  = 's3://tu-bucket/data/samplesheet.csv'
            genome       = 's3://tu-bucket/genome/genome.fa'
            gtf          = 's3://tu-bucket/genome/annotation.gtf'
            hisat2_index = 's3://tu-bucket/genome/hisat2_index'
            outdir       = 's3://tu-bucket/results'
        }

        workDir = 's3://tu-bucket/work'
    }
}
```

### Ventajas de este Enfoque Arquitectónico
1. **Desacoplamiento código/computación:** La lógica analítica (`main.nf`, módulos) permanece agnóstica a la infraestructura. Cambiar de un entorno de desarrollo local a cientos de nodos en la nube solo requiere alternar el flag de perfil (`-profile local` $\rightarrow$ `-profile aws`).
2. **Paralelismo elástico:** Muestras independientes se ejecutan de manera concurrente en instancias EC2 dedicadas (por ejemplo, `c5.2xlarge` para alineamiento y `r5.xlarge` para modelado estadístico), apagándose automáticamente al finalizar para minimizar costes.
3. **Compatibilidad con Seqera Platform:** El workflow está listo para ser monitorizado en tiempo real mediante **Nextflow Tower / Seqera Platform**, permitiendo trazabilidad de huella de carbono, consumo de memoria y observabilidad analítica.

---

## 🌍 Aplicaciones en la Práctica Biomédica

Pipelines transcriptómicos estandarizados bajo esta arquitectura son herramientas clave en:
- **Oncología de Precisión:** Identificación de sobreexpresión de oncogenes y marcadores tumorales para selección de terapias dirigidas.
- **Descubrimiento y Validación de Fármacos:** Perfilado de dianas transcripcionales en modelos preclínicos tratados con oligonucleótidos antisense o inhibidores moleculares.
- **Genómica del Neurodesarrollo:** Cuantificación precisa de perturbaciones en redes de regulación génica causadas por teratógenos o variantes genéticas raras.

---

## 📄 Licencia

Este proyecto está distribuido bajo la licencia [MIT](LICENSE), permitiendo su uso, modificación y distribución para fines académicos, biotecnológicos y profesionales.
