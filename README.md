# AASB_Portfolio_Gr06

Biblioteca Python que implementa os principais algoritmos abordados na UC.
O projeto foi desenvolvido com foco em correção algorítmica, qualidade de código, documentação completa e testes unitários, sendo facilmente importável e reutilizável por terceiros.

---

## 📦 Conteúdo do Projeto

Este portefólio inclui implementações dos seguintes tópicos:

### 1. Sequências Biológicas
- Validação de sequências (DNA, RNA e proteínas)
- Complemento e reverso-complemento
- Transcrição DNA → RNA

### 2. Alinhamento de Sequências
- Matrizes de substituição (PAM / BLOSUM)
- Needleman–Wunsch (alinhamento global)
- Smith–Waterman (alinhamento local)
- Alinhamento múltiplo progressivo
- Reconstrução dos alinhamentos e cálculo de score

### 3. Motifs e Padrões
- Procura de padrões fixos com ambiguidades
- Conversão PROSITE → expressões regulares
- Enzimas de restrição → regex + posições de corte
- PWM e PSSM:
  - Construção
  - Probabilidade de uma sequência
  - Subsequência mais provável

### 4. BLAST Simplificado
- Query map
- Identificação de hits
- Extensão dos alinhamentos
- Melhor alinhamento local

### 5. Análise Filogenética
- Construção de matriz de distâncias
- UPGMA:
  - Clustering
  - Construção da árvore filogenética


# Uso dos códigos

## Sequencias.py

```bash

#Validação de DNA
from bioinf.sequencias import validar_dna

print(validar_dna("ATGCCGTA"))   # True
print(validar_dna("atgc"))       # True
print(validar_dna("ATGX"))       # False
print(validar_dna(""))           # False

#Transcrição DNA -> RNA
from bioinf.sequencias import transcricao

dna = "ATGCCGTA"

rna = transcricao(dna)
print(rna)  # AUGCCGUA

#Complemento inverso
from bioinf.sequencias import complemento_inverso

dna = "ATGCCGTA"

comp_inv = complemento_inverso(dna)
print(comp_inv)  # TACGGCAT


```

## Alinhamento de Sequências (alinhamento.py)

```bash
# Needleman–Wunsch (Global)
from bioinf.alinhamento import needleman_wunsch

seq1 = "ATGC"
seq2 = "AGTC"

alin1, alin2, score = needleman_wunsch(seq1, seq2)

print(alin1)
print(alin2)
print("Score:", score)



# Smith–Waterman (Local)
from bioinf.alinhamento import smith_waterman

seq1 = "TGCGTAGTA"
seq2 = "CGTAGTAT"

alin1, alin2, score = smith_waterman(seq1, seq2)

print(alin1)
print(alin2)
print("Score:", score)

# Geração de consenso entre duas sequências alinhadas
from bioinf.alinhamento import gerar_consenso

alin1 = "AT-GC"
alin2 = "A-TGC"

consenso = gerar_consenso(alin1, alin2)
print(consenso)

# Alinhamento múltiplo progressivo (heurístico)
from bioinf.alinhamento import alinhamento_multiplo

sequencias = ["AC", "AT", "ATC", "CA"]

consenso_final = alinhamento_multiplo(sequencias)

print("Consenso final:", consenso_final)

```

## Blast.py

```bash
# BLAST simplificado
from bioinf.blast import blast_simplificado

query = "ATGCG"
seq_alvo = "GGATGCGTAA"

sub_q, sub_t, score, pos = blast_simplificado(query, seq_alvo)

print("Query alinhada:", sub_q)
print("Alvo alinhado:", sub_t)
print("Score:", score)
print("Posição no alvo:", pos)

# BLAST com mismatches permitidos

sub_q, sub_t, score, pos = blast_simplificado(
    query="ATGCG",
    seq_alvo="ATGAG",
    w=3,
    match=2,
    mismatch=-1
)

print(sub_q)
print(sub_t)
print(score)

```
## Filogenia.py

```bash
# Distância de Levenshtein

from bioinf.filogenia import distancia_levenshtein

print(distancia_levenshtein("ACGT", "ACGT"))   # 0
print(distancia_levenshtein("ACGT", "AGGT"))   # 1
print(distancia_levenshtein("GATTACA", "GCATGCU"))  # 4

# UPGMA
from bioinf.filogenia import upgma_simples

sequencias = ["AC", "AT", "AG", "GC"]

arvore = upgma_simples(sequencias)

print(arvore)

```
## Motifs.py

```bash
#prosite_para_regex
padrao = "C-x(2)-{GP}-[ST]"
regex = prosite_para_regex(padrao)
print(regex)

#procurar_motifs
seq = "ACCTAGCTACCT"
padrao = "C-x(2)-T"

pos = procurar_motifs(seq, padrao)
print(pos)

# criar_pwm
seqs = [
    "ATG",
    "ATC",
    "AAG",
    "ATG"
]

pwm = criar_pwm(seqs, pseudocount=1)
for i, col in enumerate(pwm):
    print(f"Posição {i}:", col)


# probabilidade_seq_pwm
seq = "ATG"
prob = probabilidade_seq_pwm(pwm, seq)
print(prob)

#subsequencia_mais_provavel
alvo = "TTATGCGATGAC"

pos, sub, prob = subsequencia_mais_provavel(pwm, alvo)
print(pos, sub, prob)
